import binascii
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import ROOT as R
import sys

from multiprocessing import Pool
import glob2
import log_loader_time
class reader:
    def __init__(self, path):
        """
        Elabora i file TL.
        """
        self.real_time=False
        self.path = path
        R.gROOT.ProcessLine('struct TreeStructTL {\
                int runNo;\
                int subRunNo;\
                int layer;\
                int gemroc;\
                int tiger;\
                int channel;\
                int tac;\
                int last_frame;\
                int tcoarse;\
                int tcoarse_10b;\
                int ecoarse;\
                int tfine;\
                int efine;\
                Double_t timestamp;\
                Double_t lcl_timestamp;\
                int delta_coarse;\
                };')
        self.time_start={}
        self.time_end={}
    def __del__(self):
        pass

    def elab_on_run(self, run,real_time):
        if real_time:
            self.real_time=real_time
            red = log_loader_time.reader(self.path, "/media/alb/Removibile/dati/")
            self.time_start,self.time_end=red.elab_on_run_dict(run)
        start_time= time.time()
        pool = Pool(processes=8)
        input_list=[]
        if not os.path.exists('/media/alb/Removibile/dati/{}'.format(run)):
            os.mkdir('/media/alb/Removibile/dati/{}'.format(run))
        for filename, (subrun, gemroc)in glob2.iglob(self.path+"/RUN_{}/SubRUN_*_GEMROC_*_TL.dat".format(run), with_matches=True):
            input_list.append((filename, gemroc, run, subrun ))
        pool.starmap(self.write_root, input_list)
        print ("All done in {:02f}".format(time.time()-start_time))

    def write_root(self, path, gemroc, run, subrun,outname="default"):
        if self.real_time:
            time_0 = int(self.time_start[subrun])
        else:
            time_0 = 0
        start_time = time.time()
        if outname=="default":
            rname = '/media/alb/Removibile/dati/{}/sub_{}_G_{}.root'.format(run,subrun,gemroc)
        else:
            rname=outname
        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStructTL()

        tree.Branch('runNo', R.AddressOf(tree_struct, 'runNo'), 'runNo/I')
        tree.Branch('subRunNo', R.AddressOf(tree_struct, 'subRunNo'), 'subRunNo/I')
        tree.Branch('layer', R.AddressOf(tree_struct, 'layer'), 'layer/I')
        tree.Branch('gemroc', R.AddressOf(tree_struct, 'gemroc'), 'gemroc/I')
        tree.Branch('tiger', R.AddressOf(tree_struct, 'tiger'), 'tiger/I')
        tree.Branch('channel', R.AddressOf(tree_struct, 'channel'), 'channel/I')
        tree.Branch('tac', R.AddressOf(tree_struct, 'tac'), 'tac/I')
        tree.Branch('last_frame', R.AddressOf(tree_struct, 'last_frame'), 'last_frame/I')
        tree.Branch('tcoarse', R.AddressOf(tree_struct, 'tcoarse'), 'tcoarse/I')
        tree.Branch('tcoarse_10b', R.AddressOf(tree_struct, 'tcoarse_10b'), 'tcoarse_10b/I')
        tree.Branch('ecoarse', R.AddressOf(tree_struct, 'ecoarse'), 'ecoarse/I')
        tree.Branch('tfine', R.AddressOf(tree_struct, 'tfine'), 'tfine/I')
        tree.Branch('efine', R.AddressOf(tree_struct, 'efine'), 'efine/I')
        tree.Branch('timestamp', R.AddressOf(tree_struct, 'timestamp'), 'timestamp/D')
        tree.Branch('lcl_timestamp', R.AddressOf(tree_struct, 'lcl_timestamp'), 'lcl_timestamp/D')

        tree.Branch('delta_coarse', R.AddressOf(tree_struct, 'delta_coarse'), 'delta_coarse/I')


        tree_struct.subRunNo = int(subrun)
        tree_struct.runNo = int(run)
        tree_struct.gemroc = int(gemroc)

        statinfo = os.stat(path)
        self.last_frame=np.zeros(8)
        self.roll_counter=np.zeros(8)

        with open(path, 'rb') as f:
            with open(rname + "missing_frames", "w") as fo:

                for i in range(0, statinfo.st_size // 8):
                    data = f.read(8)
                    if sys.version_info[0]== 2:
                        hexdata = str(binascii.hexlify(data))
                    else:
                        hexdata = str(binascii.hexlify(data), 'ascii')
                    string= "{:064b}".format(int(hexdata,16))
                    inverted=[]
                    for i in range (8,0,-1):
                        inverted.append(string[(i-1)*8:i*8])
                    string_inv="".join(inverted)
                    int_x = int(string_inv,2)
                    if (((int_x & 0xFF00000000000000) >> 59) == 0x04):  # It's a framword

                        this_framecount = ((int_x >> 15) & 0xFFFF)
                        this_tiger = ((int_x >> 56) & 0x7)
                        if self.last_frame[this_tiger]!=0 and self.last_frame[this_tiger]!=2**16-1:
                            if self.last_frame[this_tiger]!=this_framecount-1:
                                    fo.write("tiger {}, last frame {}, this frame {}, diff {}\n".format(this_tiger,self.last_frame[this_tiger],this_framecount,this_framecount-self.last_frame[this_tiger]))
                        if self.last_frame[this_tiger]!=0 and self.last_frame[this_tiger]>this_framecount:
                            self.roll_counter[this_tiger]+=1
                        self.last_frame[this_tiger]=this_framecount

                    if (((int_x & 0xFF00000000000000) >> 59) == 0x00):
                        tree_struct.tiger = (int_x>>56)&0x7
                        if self.last_frame[int(tree_struct.tiger)] != 0:
                            tree_struct.last_frame = int(self.last_frame[int(tree_struct.tiger)])
                            tree_struct.channel = (int_x>>48)&0x3F
                            tree_struct.tac = (int_x>>46)&0x3
                            tree_struct.tcoarse = (int_x >> 30)&0xFFFF
                            tree_struct.ecoarse = (int_x >> 20)&0x3FF
                            tree_struct.tfine = (int_x >> 10)&0x3FF
                            tree_struct.efine = int_x & 0x3FF
                            tree_struct.tcoarse_10b = (int_x >> 30)&0x3FF

                            if (tree_struct.last_frame%2 == 1) : # Voglio usare sempre la frame pari
                                tree_struct.lcl_timestamp = (self.roll_counter[this_tiger] * (2 ** 15) * (2 ** 16) + (tree_struct.last_frame - 1) * (2 ** 15) + tree_struct.tcoarse) * 6.25 * (10 ** (-9))
                                tree_struct.timestamp = time_0+tree_struct.lcl_timestamp
                            elif(tree_struct.last_frame%2 == 0 and tree_struct.tcoarse>2**15): #Se però arriva la pari nuova e sono sul tcoarse vecchio (>2**15) allora prendo la pari precedente
                                tree_struct.lcl_timestamp = (self.roll_counter[this_tiger] * (2 ** 15) * (2 ** 16) + (tree_struct.last_frame - 2) * (2 ** 15) + tree_struct.tcoarse) * 6.25 * (10 ** (-9))
                                tree_struct.timestamp = time_0 + tree_struct.lcl_timestamp
                            else: # Se è pari e tcoarse <2**15 prendo quella che c'è
                                tree_struct.lcl_timestamp = (self.roll_counter[this_tiger]*(2**15)*(2**16)+tree_struct.last_frame*(2**15)+tree_struct.tcoarse) * 6.25 *( 10**(-9))
                                tree_struct.timestamp = time_0+tree_struct.lcl_timestamp

                            if (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))>0:
                                tree_struct.delta_coarse = (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))
                            else:
                                tree_struct.delta_coarse = (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF)) + 1024

                            if(int(gemroc)<4):
                                tree_struct.layer = 1
                            elif(int(gemroc)<11):
                                tree_struct.layer = 2
                            else:
                                tree_struct.layer = 3
                            tree.Fill()




        rootFile.Write()
        rootFile.Close()
        print ("Sub {} G {} done in {:02f}".format(subrun, gemroc, time.time()-start_time))

    def extract_frame_in_txt(self, path, gemroc, run, subrun):
        statinfo = os.stat(path)
        self.last_frame = np.zeros(8)
        outpath="/media/alb/Removibile/dati/{}/frame_txt/subrun_{}_gemroc_{}.txt".format(run,subrun,gemroc)
        with open(path, 'rb') as f:
            with open (outpath, 'w+') as fo:
                for i in range(0, statinfo.st_size // 8):
                    data = f.read(8)
                    if sys.version_info[0] == 2:
                        hexdata = str(binascii.hexlify(data))
                    else:
                        hexdata = str(binascii.hexlify(data), 'ascii')
                    string = "{:064b}".format(int(hexdata, 16))
                    inverted = []
                    for i in range(8, 0, -1):
                        inverted.append(string[(i - 1) * 8:i * 8])
                    string_inv = "".join(inverted)
                    int_x = int(string_inv, 2)
                    if (((int_x & 0xFF00000000000000) >> 59) == 0x04):  # It's a frameword

                        this_framecount = ((int_x >> 15) & 0xFFFF)
                        this_tiger = ((int_x >> 56) & 0x7)
                        fo.write(f'tiger {this_tiger} frame {this_framecount}\n')

    def write_txt(self, path, gemroc, run, subrun,outname="default"):
        statinfo = os.stat(path)

        self.last_frame = np.zeros(8)
        if outname=="default":
            outpath="/media/alb/Removibile/dati/{}/frame_txt/subrun_{}_gemroc_{}.txt".format(run,subrun,gemroc)
        else:
            outpath=outname
        with open(path, 'rb') as f:
            with open (outpath, 'w+') as fo:
                for i in range(0, statinfo.st_size // 8):
                    data = f.read(8)
                    if sys.version_info[0] == 2:
                        hexdata = str(binascii.hexlify(data))
                    else:
                        hexdata = str(binascii.hexlify(data), 'ascii')
                    string = "{:064b}".format(int(hexdata, 16))
                    inverted = []
                    for i in range(8, 0, -1):
                        inverted.append(string[(i - 1) * 8:i * 8])
                    string_inv = "".join(inverted)
                    int_x = int(string_inv, 2)
                    if (((int_x & 0xFF00000000000000) >> 59) == 0x04):  # It's a frameword

                        this_framecount = ((int_x >> 15) & 0xFFFF)
                        this_tiger = ((int_x >> 56) & 0x7)
                        fo.write(f'tiger {this_tiger} frame {this_framecount}\n')

                    if (((int_x & 0xFF00000000000000) >> 59) == 0x04):  # It's a frameword
                        s = 'TIGER ' + '%01X: ' % ((int_x >> 56) & 0x7) + 'HB: ' + 'Framecount: %08X ' % (
                                (int_x >> 15) & 0xFFFF) + 'SEUcount: %08X\n' % (int_x & 0x7FFF)

                    if (((int_x & 0xFF00000000000000) >> 59) == 0x08):
                        s = 'TIGER ' + '%01X: ' % ((int_x >> 56) & 0x7) + 'CW: ' + 'ChID: %02X ' % (
                                (int_x >> 24) & 0x3F) + ' CounterWord: %016X\n' % (int_x & 0x00FFFFFF)
                    if (((int_x & 0xFF00000000000000) >> 59) == 0x00):
                        s = 'TIGER ' + '%01X: ' % ((int_x >> 56) & 0x7) + 'EW: ' + 'ChID: %02X ' % (
                                (int_x >> 48) & 0x3F) + 'tacID: %01X ' % ((int_x >> 46) & 0x3) + 'Tcoarse: %04X ' % (
                                    (int_x >> 30) & 0xFFFF) + 'Ecoarse: %03X ' % (
                                    (int_x >> 20) & 0x3FF) + 'Tfine: %03X ' % ((int_x >> 10) & 0x3FF) + 'Efine: %03X \n' % (
                                    int_x & 0x3FF)
                    fo.write(s)




if __name__ == "__main__":
    runner = reader("/home/alb/srv_lab_raw/")
    # runner = reader("/media/alb/Removibile/dati_raw/")

    runner.elab_on_run(402,True)
    #runner.extract_frame_in_txt("/media/alb/space/TIGER_scriptsV3/data_folder/RUN_398/SubRUN_3_GEMROC_0_TL.dat", 0, 398, 3)
