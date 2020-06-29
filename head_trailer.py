import binascii
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import ROOT as R
import sys
import subprocess
from multiprocessing import Pool
import glob2

def maxDiff(a):
    vmin = a[0]
    dmax = 0
    for i in range(len(a)):
        if (a[i] < vmin):
            vmin = a[i]
        elif (a[i] - vmin > dmax):
            dmax = a[i] - vmin
    return dmax
class reader:
    def __init__(self, inpath, outpath):
        self.path = inpath
        self.outpath = outpath
        R.gROOT.ProcessLine('struct TreeStruct {\
                int runNo;\
                int subRunNo;\
                int gemroc;\
                int count;\
                int l1_ts;\
                int hit_count;\
                int last_l1_ts_dif;\
                int l1_frame;\
                int tiger_id;\
                int UDP_num;\
                int count_trailer;\
                int ch;\
                int last_count_from_ch;\
                bool top_L1_chk_error;\
                bool header_misalignment_error;\
                bool FIFO_FULL_error;\
                bool daq_pll_unlocked;\
                bool global_rx_error;\
                bool XCVR_rx_alignment_error;\
                int no_trailer_bug;\
        };')
    def __del__(self):
        pass

    def run_on_run(self, run):
        start_time= time.time()
        pool = Pool(processes=8)
        input_list=[]
        if not os.path.exists('{}/{}'.format(self.outpath, run)):
            print ("Making folder for decoded header and trailer info")
            os.mkdir('{}/{}'.format(self.outpath, run))
        for filename, (subrun, gemroc) in glob2.iglob(self.path + "/RUN_{}/SubRUN_*_GEMROC_*_TM.dat".format(run), with_matches=True):
            input_list.append((filename, gemroc, run, subrun ))
        print ("{} raw data file found".format(len(input_list)))
        pool.starmap(self._write_root, input_list)
        print ("All done in {:02f}".format(time.time()-start_time))
        # subprocess.call(['/bin/bash', '-c', "hadd -f {0}/{1}.root {0}/{1}/sub_*_G_*".format(self.outpath,run)])
        # for filename in glob2.glob("{}/{}/sub_*_G_*".format(self.outpath,run)):
        #     os.remove(filename)
        # print ("Merging done")


    def _write_root(self, path, gemroc, run, subrun):
        start_time = time.time()

        rname = '{}/{}/sub_{}_G_{}.root'.format(self.outpath,run,subrun,gemroc)


        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStruct()
        tree.Branch('runNo', R.AddressOf(tree_struct, 'runNo'), 'runNo/I')
        tree.Branch('subRunNo', R.AddressOf(tree_struct, 'subRunNo'), 'subRunNo/I')
        tree.Branch('gemroc', R.AddressOf(tree_struct, 'gemroc'), 'gemroc/I')
        tree.Branch('count', R.AddressOf(tree_struct, 'count'), 'count/I')
        tree.Branch('l1_ts', R.AddressOf(tree_struct, 'l1_ts'), 'l1_ts/I')
        tree.Branch('hit_count', R.AddressOf(tree_struct, 'hit_count'), 'hit_count/I')
        tree.Branch('last_l1_ts_dif', R.AddressOf(tree_struct, 'last_l1_ts_dif'), 'last_l1_ts_dif/I')
        tree.Branch('l1_frame', R.AddressOf(tree_struct, 'l1_frame'), 'l1_frame/I')
        tree.Branch('tiger_id', R.AddressOf(tree_struct, 'tiger_id'), 'tiger_id/I')
        tree.Branch('count_trailer', R.AddressOf(tree_struct, 'count_trailer'), 'count_trailer/I')
        tree.Branch('ch', R.AddressOf(tree_struct, 'ch'), 'ch/I')
        tree.Branch('last_count_from_ch', R.AddressOf(tree_struct, 'last_count_from_ch'), 'last_count_from_ch/I')
        tree.Branch('UDP_num', R.AddressOf(tree_struct, 'UDP_num'), 'UDP_num/I')
        tree.Branch('top_L1_chk_error', R.AddressOf(tree_struct, 'top_L1_chk_error'), 'top_L1_chk_error/b')
        tree.Branch('header_misalignment_error', R.AddressOf(tree_struct, 'header_misalignment_error'), 'header_misalignment_error/b')
        tree.Branch('FIFO_FULL_error', R.AddressOf(tree_struct, 'FIFO_FULL_error'), 'FIFO_FULL_error/b')
        tree.Branch('daq_pll_unlocked', R.AddressOf(tree_struct, 'daq_pll_unlocked'), 'daq_pll_unlocked/b')
        tree.Branch('global_rx_error', R.AddressOf(tree_struct, 'global_rx_error'), 'global_rx_error/b')
        tree.Branch('XCVR_rx_alignment_error', R.AddressOf(tree_struct, 'XCVR_rx_alignment_error'), 'XCVR_rx_alignment_error/b')
        tree.Branch('no_trailer_bug', R.AddressOf(tree_struct, 'no_trailer_bug'), 'no_trailer_bug/I')


        tree_struct.subRunNo = int(subrun)
        tree_struct.runNo = int(run)
        tree_struct.gemroc = int(gemroc)


        statinfo = os.stat(path)
        self.last_frame=np.zeros(8)

        last_l1_ts=-1

        done_header = False
        done_trailer = False
        done_UDP = False
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

                if (((int_x & 0xE000000000000000) >> 61) == 0x6): #HEADER
                    LOCAL_L1_COUNT_31_6 = int_x >> 32 & 0x3FFFFFF
                    LOCAL_L1_COUNT_5_0 = int_x >> 24 & 0x3F
                    LOCAL_L1_COUNT = (LOCAL_L1_COUNT_31_6 << 6) + LOCAL_L1_COUNT_5_0
                    LOCAL_L1_TIMESTAMP = int_x & 0xFFFF
                    HITCOUNT = (int_x >> 16) & 0xFF
                    tree_struct.count = LOCAL_L1_COUNT
                    tree_struct.l1_ts = LOCAL_L1_TIMESTAMP
                    tree_struct.hit_count = HITCOUNT

                    tree_struct.top_L1_chk_error = (int_x >> 58)& 0x1
                    tree_struct.header_misalignment_error = (int_x >> 59)& 0x1
                    tree_struct.FIFO_FULL_error = (int_x >> 60)& 0x1

                    if last_l1_ts!=-1:
                        if (((int_x & 0xFFFF) - last_l1_ts) > 0):
                            L1_TS_abs_diff = LOCAL_L1_TIMESTAMP - last_l1_ts
                        else:
                            L1_TS_abs_diff = 65536 + LOCAL_L1_TIMESTAMP - last_l1_ts
                        tree_struct.last_l1_ts_dif = L1_TS_abs_diff
                    last_l1_ts = LOCAL_L1_TIMESTAMP
                    done_header= True
                if (((int_x & 0xE000000000000000) >> 61) == 0x7): #TRAILER
                    tree_struct.l1_frame = ((int_x >> 37) & 0xFFFFFF)
                    tree_struct.tiger_id = ((int_x >> 27) & 0x7)
                    tree_struct.count_trailer = ((int_x >> 24) & 0x7)
                    tree_struct.ch = ((int_x >> 18) & 0x3F)
                    tree_struct.last_count_from_ch = (int_x & 0x3FFFF)
                    done_trailer= True

                if (((int_x & 0xF000000000000000) >> 60) == 0x4): #UDP trailer
                    tree_struct.UDP_num = ((int_x >> 32) & 0xFFFFF) + ((int_x >> 0) & 0xFFFFFFF)
                    tree_struct.daq_pll_unlocked = (int_x >> 57)& 0x1
                    tree_struct.global_rx_error = (int_x >> 58)& 0x1
                    tree_struct.XCVR_rx_alignment_error = (int_x >> 59)& 0x1
                    if not done_header:
                        print ("missing header!")
                        tree_struct.no_trailer_bug=3
                    if not done_trailer:
                        tree_struct.no_trailer_bug=1
                    else:
                        tree_struct.no_trailer_bug=0

                    tree.Fill()
                    done_header = False
                    done_trailer = False



        rootFile.Write()
        rootFile.Close()
        # print ("Sub {} G {} done in {:02f}".format(subrun, gemroc, time.time()-start_time))

class analyzer:
    def __init__(self, path, outpath):
        self.inpath = path
        self.outpath = outpath
        self.run=0
        R.gROOT.ProcessLine('struct TreeStruct_2 {\
                int runNo;\
                int subRunNo;\
                int count;\
                int l1_ts_diff;\
        };')
    def write_root_time_dif(self):
        rname = '{}/{}_time_diff.root'.format(self.outpath, self.run)

        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStruct_2()
        tree.Branch('runNo', R.AddressOf(tree_struct, 'runNo'), 'runNo/I')
        tree.Branch('subRunNo', R.AddressOf(tree_struct, 'subRunNo'), 'subRunNo/I')
        tree.Branch('count', R.AddressOf(tree_struct, 'count'), 'count/I')
        tree.Branch('l1_ts_diff', R.AddressOf(tree_struct, 'l1_ts_diff'), 'l1_ts_diff/I')
        f = R.TFile.Open('{}/{}.root'.format(self.outpath, self.run))
        tree_struct.runNo = int(self.run)

        trigger_dict = {}
        for entryNum in range(0, f.tree.GetEntries()):
            f.tree.GetEntry(entryNum)
            gemroc = int(getattr(f.tree, "gemroc"))
            subRunNo = int(getattr(f.tree, "subRunNo"))
            count = int(getattr(f.tree, "count"))
            l1_ts = int(getattr(f.tree, "l1_ts"))
            if subRunNo not in trigger_dict:
                trigger_dict[subRunNo] = {}
            if count not in trigger_dict[subRunNo]:
                trigger_dict[subRunNo][count] = {}
            trigger_dict[subRunNo][count][gemroc] = l1_ts
        for subrun in trigger_dict:
            for count in trigger_dict[subrun]:
                listina = []
                for GEMROC in trigger_dict[subrun][count]:
                    listina.append(trigger_dict[subrun][count][GEMROC])
                tree_struct.subRunNo = subrun
                tree_struct.count = count
                tree_struct.l1_ts_diff = maxDiff(listina)
                tree.Fill()
        rootFile.Write()
        rootFile.Close()

    def build_trig_dict(self, chain):
        trigger_dict={}
        for  entryNum  in  range(0,chain.GetEntries ()):
            chain.GetEntry(entryNum)
            count=getattr(chain,"count")
            gemroc=getattr(chain,"gemroc")
            if gemroc not in trigger_dict.keys():
                trigger_dict[gemroc]= []
            if count not in trigger_dict[gemroc]:
                trigger_dict[gemroc].append(count)
        return trigger_dict



    def load_subrun(self, sub):
        ch = R.TChain("tree")
        (ch.Add("{}/{}/sub_{}_G*".format(self.inpath, self.run, sub)))
        return ch




    def find_missing(self, trigger_dict, sub):
        for gemroc in trigger_dict:
            missing=(missing_elements(trigger_dict[gemroc],0, len(trigger_dict[gemroc]) - 1))
            for elem in missing:
                with open (self.outpath+"/skipped_packets_run_{}".format(self.run),"a+") as fo:
                    fo.write("subrun {} : gemroc {} missing {}\n".format(sub, gemroc, elem))
                print ("subrun {} : gemroc {} missing {}".format(sub, gemroc, elem))



    def run_on_run(self, run):
        self.run=run
        if os.path.exists(self.outpath+"/skipped_packets_run_{}".format(self.run)):
            os.remove(self.outpath+"/skipped_packets_run_{}".format(self.run))
        sub_list=[]
        pool = Pool(processes=8)
        start_time=time.time()
        print (" {}/{}/sub_*_G*".format(self.inpath, run))
        for filename,(sub,roc) in glob2.iglob("{}/{}/sub_*_G*".format(self.inpath, run), with_matches=True):
            if sub not in sub_list:
                sub_list.append(sub)
        print ("{} subrun found ".format(len (sub_list)))
        pool.map(self.process_data, sub_list)
        print ("All done in {:02f}".format(time.time()-start_time))

    def process_data (self, sub):
        ch = self.load_subrun(sub)
        trigger_dict = self.build_trig_dict(ch)
        self.find_missing(trigger_dict, sub)

def missing_elements(L, start, end):
    if end - start <= 1:
        if L[end] - L[start] > 1:
            yield from range(L[start] + 1, L[end])
        return

    index = start + (end - start) // 2

    # is the lower half consecutive?
    consecutive_low = L[index] == L[start] + (index - start)
    if not consecutive_low:
        yield from missing_elements(L, start, index)

    # is the upper part consecutive?
    consecutive_high = L[index] == L[end] - (end - index)
    if not consecutive_high:
        yield from missing_elements(L, index, end)


if __name__ == "__main__":
    # reader_1=reader("/home/alb/srv_lab_raw", "/home/alb/Desktop/elaborazioni_e_dati/analisi_run/trigger_decode")
    # reader_1.run_on_run(405)
    analyzer_1=analyzer("/home/alb/Desktop/elaborazioni_e_dati/analisi_run/trigger_decode","/home/alb/Desktop/elaborazioni_e_dati/analisi_run/out")
    analyzer_1.run_on_run(405)