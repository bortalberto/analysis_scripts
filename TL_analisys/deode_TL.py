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

class reader:
    def __init__(self, path):
        self.path = path
        R.gROOT.ProcessLine('struct TreeStruct {\
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
                float timestamp;\
                int delta_coarse;\
                };')
    def __del__(self):
        pass

    def elab_on_run(self, run):
        start_time= time.time()
        pool = Pool(processes=8)
        input_list=[]
        if not os.path.exists('/home/alb/Desktop/elaborazioni_e_dati/analisi_run/TL_analisys/data_out/{}'.format(run)):
            os.mkdir('/home/alb/Desktop/elaborazioni_e_dati/analisi_run/TL_analisys/data_out/{}'.format(run))
        for filename, (subrun, gemroc)in glob2.iglob(self.path+"/RUN_{}/SubRUN_*_GEMROC_*_TL.dat".format(run), with_matches=True):
            input_list.append((filename, gemroc, run, subrun ))
        pool.starmap(self.write_root, input_list)
        print ("All done in {:02f}".format(time.time()-start_time))

    def write_root(self, path, gemroc, run, subrun):
        start_time = time.time()

        rname = '/home/alb/Desktop/elaborazioni_e_dati/analisi_run/TL_analisys/data_out/{}/sub_{}_G_{}.root'.format(run,subrun,gemroc)


        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStruct()
        tree_struct.subRunNo = int(subrun)
        tree_struct.runNo = int(run)
        tree_struct.gemroc = int(gemroc)

        for key in R.TreeStruct.__dict__.keys():
            if '__' not in key:
                formstring = '/F'
                if isinstance(tree_struct.__getattribute__(key), int):
                    formstring = '/I'
                tree.Branch(key, R.AddressOf(tree_struct, key), key + formstring)


        statinfo = os.stat(path)
        self.last_frame=np.zeros(8)

        with open(path, 'rb') as f:
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
                        if tree_struct.last_frame//2==0 and tree_struct.tcoarse > 2**15:
                            tree_struct.timestamp = ((tree_struct.last_frame-1)*2**15+tree_struct.tcoarse) *6.25 * 10**(-9)
                        else:
                            tree_struct.timestamp = (tree_struct.last_frame*2**15+tree_struct.tcoarse) * 6.25 * 10**(-9)
                        if (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))>0:
                            tree_struct.delta_coarse = (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))
                        else:
                            tree_struct.delta_coarse = (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF)) + 1024
                        temp_ecoarse = tree_struct.ecoarse
                        tree_struct.charge_SH = int_x & 0x3FF

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





if __name__ == "__main__":
    runner = reader("/media/alb/space/TIGER_scriptsV3/data_folder")
    runner.elab_on_run(398)
    #runner.write_root("/media/alb/space/TIGER_scriptsV3/data_folder/RUN_398/SubRUN_3_GEMROC_0_TL.dat", 0, 398, 3)
