import glob2
from multiprocessing import Pool
import os
import time
import numpy as np
import ROOT as R
import subprocess
import time
import sys
from calendar import timegm
import pickle


class reader:
    def __init__(self, path, elabpath):
        self.path = path
        self.elabpath=elabpath
        R.gROOT.ProcessLine('struct TreeStruct {\
                int runNo;\
                int subRunNo;\
                int triggers;\
                int start_time;\
                int end_time;\
        };')

    def elab_on_run(self, run):
        rname= self.elabpath+"{}".format(run)
        rname = rname + '/time_run_{}.root'.format(run)


        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStruct()

        tree.Branch('runNo', R.AddressOf(tree_struct, 'runNo'), 'runNo/I')
        tree.Branch('subRunNo', R.AddressOf(tree_struct, 'subRunNo'), 'subRunNo/I')
        tree.Branch('triggers', R.AddressOf(tree_struct, 'triggers'), 'triggers/I')
        tree.Branch('start_time', R.AddressOf(tree_struct, 'start_time'), 'start_time/I')
        tree.Branch('end_time', R.AddressOf(tree_struct, 'end_time'), 'end_time/I')
        tree_struct.runNo = int(run)

        for filename, (subrun,)in glob2.iglob(self.path+"/RUN_{}/ACQ_log_*".format(run), with_matches=True):
            tree_struct.subRunNo = int(subrun)
            lenght_list = []
            with open(filename, "r") as flog:
                lines = flog.readlines()
            for line in lines:
                if "packets" in line:
                    lenght_list.append (int (line.split("TM.dat, total packets= ")[1]))
            if len(lines)!=0:
                if "--" in lines[0]:
                    utc_time = time.strptime(lines[0].split("--")[0], "%a %b  %d %H:%M:%S %Y  ")
                    start_time = timegm(utc_time)
                else:
                    start_time = -1
                if "--" in lines[-1]:
                    utc_time = time.strptime(lines[-1].split("--")[0], "%a %b  %d %H:%M:%S %Y  ")
                    end_time= timegm(utc_time)
                else:
                    end_time = start_time
                tree_struct.start_time=start_time
                tree_struct.end_time=end_time
            else:
                tree_struct.start_time=-1
                tree_struct.end_time=-1
            lenght_list=np.asarray(lenght_list)

            if len(lenght_list)>0 and (np.count_nonzero(lenght_list)): #interrotto da timeout, con dati
                tree_struct.triggers = int (np.min(lenght_list[np.nonzero(lenght_list)]))
            else:
                tree_struct.triggers=-1

            tree.Fill()
        rootFile.Write()
        rootFile.Close()

    def elab_on_run_dict(self, run):
        time_dict_start={}
        time_dict_end={}

        for filename, (subrun,)in glob2.iglob(self.path+"/RUN_{}/ACQ_log_*".format(run), with_matches=True):
            lenght_list = []
            with open(filename, "r") as flog:
                lines = flog.readlines()
            if len(lines)!=0:
                if "--" in lines[0]:
                    utc_time = time.strptime(lines[0].split("--")[0], "%a %b  %d %H:%M:%S %Y  ")
                    start_time = timegm(utc_time)
                else:
                    start_time = -1
                if "--" in lines[-1]:
                    utc_time = time.strptime(lines[-1].split("--")[0], "%a %b  %d %H:%M:%S %Y  ")
                    end_time= timegm(utc_time)
                else:
                    end_time = start_time
                time_dict_start[subrun]=start_time
                time_dict_end[subrun]=end_time

            else:
                time_dict_start[subrun]=-1
                time_dict_end[subrun]=-1

        return time_dict_start, time_dict_end

    def save_dict(self, time_dict,run):
        with open(self.elabpath+"/"+run, "wb+") as savefile:
            pickle.dump(time_dict, savefile)
if __name__ == "__main__":
        if len (sys.argv)!=2:
            print ("Time extraction")
        else:
            run = sys.argv[1]
            # runner = reader(os.environ["data"]+"/raw_dat", os.environ["data"]+"/raw_dat" )
            runner = reader("/home/alb/srv_lab_raw/", "/media/alb/Removibile/time_save/" )
            # runner.elab_on_run(run)
            time_dict=runner.elab_on_run_dict(run)
            runner.save_dict(time_dict,run)
            print("Time extraction done")





