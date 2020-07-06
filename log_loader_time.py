import glob2
from multiprocessing import Pool
import os
import time
import numpy as np
import ROOT as R
import subprocess

import sys
class reader:
    def __init__(self, path, elabpath, outpath):
        self.path = path
        self.elabpath=elabpath
        self.outpath=outpath
        R.gROOT.ProcessLine('struct TreeStruct {\
                int runNo;\
                int subRunNo;\
                int triggers;\
                int error_end_FEB;\
                int standard_end;\
                int timeout_end;\
                int DAQ_crash;\
                int end;\
        };')

    def elab_on_run(self, run):
        start_time= time.time()
        pool = Pool(processes=8)
        input_list=[]

        for filename, (subrun)in glob2.iglob(self.path+"/RUN_{}/ACQ_log_*".format(run), with_matches=True):
            input_list.append((filename, run, subrun[0]))
        rname= self.elabpath +"/"+ run + "/end_of_run"
        if not os.path.exists(rname):
            os.mkdir(rname)
        pool.starmap(self.extract_stop_reason, input_list)
        print ("All done in {:02f}".format(time.time()-start_time))
        subprocess.call(['/bin/bash', '-i', '-c', "hadd -f {0}/run_{1} {0}/run_{1}sub_*".format(self.elabpath +"/"+ run + "/end_of_run", run)])
        for filename in glob2.glob("{}/run_{}sub*".format(self.elabpath +"/"+ run + "/end_of_run", run)):
            os.remove(filename)

    def extract_stop_reason(self, filename, run, subrun):
        rname= self.elabpath +"/"+ run + "/end_of_run"
        rname = rname + '/run_{}sub_{}.root'.format(run, subrun)


        rootFile = R.TFile(rname, 'recreate')
        tree = R.TTree('tree', '')
        tree_struct = R.TreeStruct()

        tree.Branch('runNo', R.AddressOf(tree_struct, 'runNo'), 'runNo/I')
        tree.Branch('subRunNo', R.AddressOf(tree_struct, 'subRunNo'), 'subRunNo/I')
        tree.Branch('triggers', R.AddressOf(tree_struct, 'triggers'), 'triggers/I')
        tree.Branch('error_end_FEB', R.AddressOf(tree_struct, 'error_end_FEB'), 'error_end_FEB/I')
        tree.Branch('standard_end', R.AddressOf(tree_struct, 'standard_end'), 'standard_end/I')
        tree.Branch('timeout_end', R.AddressOf(tree_struct, 'timeout_end'), 'timeout_end/I')
        tree.Branch('DAQ_crash', R.AddressOf(tree_struct, 'DAQ_crash'), 'DAQ_crash/I')
        tree.Branch('end', R.AddressOf(tree_struct, 'end'), 'end/I')

        tree_struct.subRunNo = int(subrun)
        tree_struct.runNo = int(run)
        lenght_list = []
        errors_list = []
        timeout = False
        dct_list=[(8,4),(8,5),(12,4),(12,5),(12,6),(12,7)]
        with open(filename, "r") as flog:
            lines = flog.readlines()
        for line in lines:
            if "packets" in line:
                lenght_list.append (int (line.split("TM.dat, total packets= ")[1]))
            if " 16777215 8/10 bit errors" in line:
                GEMROC = int(line.split("GEMROC")[1].split(" ")[1])
                TIGER = int(line.split("TIGER")[1].split(" ")[1])
                errors_list.append((GEMROC, TIGER))
            if "out" in line:
                timeout = True
        for dontcare in dct_list:
            if dontcare in errors_list: errors_list.remove(dontcare)
        lenght_list=np.asarray(lenght_list)
        # if len (lenght_list)>0:
        #     print ("Run {} subrun {} triggers: {}".format (run, int(subrun), np.average(lenght_list)))
        # else:
        #     print ("Run {} subrun {} never started or DAQ crashed".format (run, int(subrun)))
        if timeout:
            if len(lenght_list)>0 and (np.count_nonzero(lenght_list)): #interrotto da timeout, con dati
                # print ("Run {} subrun {} ended for timeout triggers: {}".format (run, int(subrun), np.min(lenght_list)))
                # tree_struct.triggers = int (np.min( list(i for i in lenght_list if i>0)))
                tree_struct.triggers = int (np.min(lenght_list[np.nonzero(lenght_list)]))

                tree_struct.standard_end = 0
                tree_struct.timeout_end = 1
                tree_struct.error_end_FEB = -1
                tree_struct.DAQ_crash = 0
                tree_struct.end = 1


        elif len (errors_list)>0: #interrotto da errori
            # print ("Run {} subrun {} ended for errors on G{} T{} triggers: {}".format (run, int(subrun),errors_list[0][0],errors_list[0][1], np.min(lenght_list)))
            tree_struct.triggers = int (np.min(lenght_list[np.nonzero(lenght_list)]))
            tree_struct.standard_end = 0
            tree_struct.timeout_end = 0
            tree_struct.error_end_FEB = errors_list[0][0]*4+errors_list[0][1]//2
            tree_struct.DAQ_crash = 0
            tree_struct.end = 2

        elif len (lenght_list)>0:#ok
            # print ("Run {} subrun {} triggers: {}".format (run, int(subrun), np.min(lenght_list)))
            tree_struct.triggers = int (np.min(lenght_list[np.nonzero(lenght_list)]))
            tree_struct.standard_end = 1
            tree_struct.timeout_end = 0
            tree_struct.error_end_FEB = -1
            tree_struct.DAQ_crash = 0
            tree_struct.end = 0
        else: #qualcosa di strano
            # print ("Run {} subrun {} never started or DAQ crashed".format (run, int(subrun)))
            tree_struct.triggers = -1
            tree_struct.standard_end = 0
            tree_struct.timeout_end = 0
            tree_struct.error_end_FEB = -1
            tree_struct.DAQ_crash = 1
            tree_struct.end = 3

        tree.Fill()
        rootFile.Write()
        rootFile.Close()
    def produce_plot(self, run):
        R.gStyle.SetOptStat(0)
        R.gROOT.SetBatch(1)

        f = R.TFile.Open(self.elabpath +"/"+ run + "/end_of_run"+"/run_{}".format(run))
        c = R.TCanvas("c")

        leg = R.TLegend(0.58, 0.7, 0.9, 0.9)

        f.tree.Draw("end")
        leg.AddEntry("", "0 = normally ended", "")
        leg.AddEntry("", "1 = timeout", "")
        leg.AddEntry("", "2 = 8b/10b errors", "")
        leg.AddEntry("", "3 = crash", "")
        leg.Draw()
        c.Draw()
        c.SaveAs(self.outpath + "/" + run + "/run_{}_sub_end.png".format(run))


if __name__ == "__main__":
        if len (sys.argv)!=2:
            print ("End of run analysis: Wrong input")
        else:
            run = sys.argv[1]
            runner = reader(os.environ["data"]+"/raw_dat", os.environ["QAQC_out"], os.environ["QAQC_show"])
            runner.elab_on_run(run)
            runner.produce_plot(run)
            print("End of run analysis done")





