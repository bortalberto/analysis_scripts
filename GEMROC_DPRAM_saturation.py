#!/usr/bin/env python
# coding: utf-8

### Riempimento della DPRAM


import ROOT as R
import glob2
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from array import array
import sys

class DPRAM_filling_check():

    def __init__(self, run):
        self.run = run
        self.tcoarse_dict = {}
        self.l1_trigger_dict = {}
        self.l1_trigger_extremes = {}
        self.length_dict = {}
        self.sub_run = 0
        self.risky = []
        self.sure_sat = []
        self.wierd = []
        self.total_trig = 0
        self.total_count_risky = 0
        self.total_count_sure_sat = 0
        self.total_count_wierd = 0

    def elab(self):
        """
        Runs the other functions on the whole run
        :return:
        """
        print("/home/alb/srv_lab/{}/Sub_RUN_ana_*".format(self.run))
        self.create_tree()
        for filename, (subrun, )in glob2.iglob("/home/alb/srv_lab/{}/Sub_RUN_ana_*.root".format(self.run), with_matches=True):
            print (filename)
            f = self.load_file(subrun)
            self.sub_run = subrun
            self.build_tcoarse_dict_and_trigger_dict(f)
            self.build_trigger_extremes()
            self.build_length_dict()
            self.analyze_data()
            self.refresh()
            self.results()
        self.close_tree()

    def elab_on_sub(self,sub):
        """
        Runs the other functions on the whole run
        :return:
        """
        print("/home/alb/srv_lab/{}/Sub_RUN_ana_*".format(self.run))
        self.create_tree()
        f = self.load_file(sub)
        self.sub_run = sub
        self.build_tcoarse_dict_and_trigger_dict(f)
        self.build_trigger_extremes()
        self.build_length_dict()
        self.analyze_data()
        self.refresh()
        self.close_tree()
        self.results()

    def load_file(self, subrun):
        filename="/home/alb/srv_lab/{}/Sub_RUN_ana_{}.root".format(self.run,subrun)
        print(f"Opening {filename}")
        return R.TFile.Open (filename)

    def create_tree(self):
        R.gROOT.ProcessLine(
        "struct tree_struct_t {\
            Int_t           run;\
            Int_t           subrun;\
            Int_t           count;\
            Int_t           gemroc;\
            Int_t           feb;\
            Int_t           bucket;\
            Int_t           n_hit;\
        };" );

        self.tree_struct = R.tree_struct_t()

        self.fo = R.TFile('DPRAM_study_run_{}.root'.format(self.run), 'recreate')
        self.to = R.TTree('tree', 'Bucket filling')
        self.to.Branch('run', R.AddressOf(self.tree_struct, 'run'),"run/I")
        self.to.Branch('subrun', R.AddressOf(self.tree_struct, 'subrun'),"subrun/I")
        self.to.Branch('count', R.AddressOf(self.tree_struct, 'count'),"count/I")
        self.to.Branch('gemroc', R.AddressOf(self.tree_struct, 'gemroc'),"gemroc/I")
        self.to.Branch('feb', R.AddressOf(self.tree_struct, 'feb'),"feb/I")
        self.to.Branch('bucket', R.AddressOf(self.tree_struct, 'bucket'),"bucket/I")
        self.to.Branch('n_hit', R.AddressOf(self.tree_struct, 'n_hit'),"n_hit/I")
        pass

    def close_tree(self):
        self.fo.Write()
        self.fo.Close()

    def refresh(self):
        """
        Clears for the next subrun
        :return:
        """
        self.tcoarse_dict = {}
        self.l1_trigger_dict = {}
        self.l1_trigger_extremes = {}
        self.length_dict = {}

    def build_tcoarse_dict_and_trigger_dict(self, fil):
        """
        Build the tcoarse dict:
            struct : [count][gemroc][FEB_group][page]
            Page = Tcoarse[8:11], ma in python [4:8]
            
        Trigger_dict:
            coarse value of the trigger
        :param fil: 
        :return: 
        """
        for entryNum in range(0, fil.tree.GetEntries()):
            fil.tree.GetEntry(entryNum)
            tcoarse = int(getattr(fil.tree, "tcoarse"))
            l1ts_min_tcoarse = int(getattr(fil.tree, "l1ts_min_tcoarse"))

            page = (f'{tcoarse:016b}'[4:8])
            count = getattr(fil.tree, "count")
            gemroc = getattr(fil.tree, "gemroc")
            tiger = getattr(fil.tree, "tiger")
            ecoarse = getattr(fil.tree, "ecoarse")
            if ecoarse == 0 and tcoarse == 0:
                # Data corrupted by the 0s instead of counter word bug
                # print ("Corrupted data")
                pass
            else:
                if tiger < 2:
                    FEB_group = 0
                elif tiger < 4:
                    FEB_group = 1
                elif tiger < 6:
                    FEB_group = 2
                elif tiger < 8:
                    FEB_group = 3

                if count not in self.tcoarse_dict.keys():
                    self.tcoarse_dict[count] = {}
                if gemroc not in self.tcoarse_dict[count].keys():
                    self.tcoarse_dict[count][gemroc] = {}
                if FEB_group not in self.tcoarse_dict[count][gemroc].keys():
                    self.tcoarse_dict[count][gemroc][FEB_group] = {}
                if page not in self.tcoarse_dict[count][gemroc][FEB_group].keys():
                    self.tcoarse_dict[count][gemroc][FEB_group][page] = 0
                self.tcoarse_dict[count][gemroc][FEB_group][page] += 1

                if count not in self.l1_trigger_dict.keys():
                    self.l1_trigger_dict[count] = []
                self.l1_trigger_dict[count].append(l1ts_min_tcoarse + tcoarse)

    def build_trigger_extremes(self):
        """
        Finds the beginning and the end of the trigger windows, taking the mode of the triggers for the same L1 count
        """
        for count in self.l1_trigger_dict:
            #     print (f'Count: {count}, Std:{np.std(l1_trigger_dict[count])}\n {l1_trigger_dict[count]}')
            self.l1_trigger_dict[count] = (stats.mode(self.l1_trigger_dict[count])).mode[0]
            if self.l1_trigger_dict[count] < 1600:
                self.l1_trigger_extremes[count] = (int(f'{self.l1_trigger_dict[count] + 4096 - 1568:016b}'[4:16], 2), int(f'{self.l1_trigger_dict[count] + 4096 - 1298:016b}'[4:16], 2))
            else:
                self.l1_trigger_extremes[count] = (int(f'{self.l1_trigger_dict[count] - 1568:016b}'[4:16], 2), int(f'{self.l1_trigger_dict[count] - 1298:016b}'[4:16], 2))

            # self.l1_trigger_extremes[count] = (int(f'{self.l1_trigger_dict[count]:016b}'[4:16], 2) - 1568, int(f'{self.l1_trigger_dict[count]:016b}'[4:16], 2) - 1298)

    def build_length_dict(self):
        """
        Build a dict containing the % of the word taken into the trigger
        :param self:
        :return:
        """
        for count in self.l1_trigger_dict:
            self.length_dict[count] = {}
            #         print (count)
            len3 = None
            beginning = (int(f'{self.l1_trigger_extremes[count][0]:016b}'[4:16], 2))
            end = (int(f'{self.l1_trigger_extremes[count][1]:016b}'[4:16], 2))
            first_page = (int(f'{self.l1_trigger_extremes[count][0]:016b}'[4:8], 2)) * 256
            last_page = (int(f'{self.l1_trigger_extremes[count][1]:016b}'[4:8], 2)) * 256
            if (int(f'{self.l1_trigger_extremes[count][1]:016b}'[4:8], 2) - int(f'{self.l1_trigger_extremes[count][0]:016b}'[4:8], 2)) in (2, -14):
                len1 = 256 - (beginning - first_page)
                len2 = 256
                len3 = end - last_page
            elif (last_page - beginning > 0):
                len1 = last_page - beginning
                len2 = end - last_page
            elif (last_page - beginning < 0):
                len1 = last_page + 4096 - beginning
                len2 = end - last_page
            #         print (len1,len2,len3,len1+len2+len3)
            if len3 == None:
                self.length_dict[count][(f'{self.l1_trigger_extremes[count][0]:016b}'[4:8])] = len1 / 256
                self.length_dict[count][(f'{self.l1_trigger_extremes[count][1]:016b}'[4:8])] = len2 / 256
            else:
                self.length_dict[count][(f'{self.l1_trigger_extremes[count][0]:016b}'[4:8])] = len1 / 256
                self.length_dict[count][(f'{self.l1_trigger_extremes[count][0] + 256:016b}'[4:8])] = len2 / 256
                self.length_dict[count][(f'{self.l1_trigger_extremes[count][1]:016b}'[4:8])] = len3 / 256

    def analyze_data(self):
        risky = []
        sure_sat = []
        wierd = []
        count_min = min(self.tcoarse_dict.keys())
        for count in self.tcoarse_dict.keys():
            if count > count_min + 3:  # Scarto i primi 3 triggers
                self.total_trig += 1
                for gemroc in self.tcoarse_dict[count].keys():
                    if gemroc in range(0,14):
                        for FEB_group in self.tcoarse_dict[count][gemroc].keys():
                            for page in self.tcoarse_dict[count][gemroc][FEB_group].keys():
                                try:
                                    if (self.tcoarse_dict[count][gemroc][FEB_group][page] >= self.length_dict[count][page] * 32) :
                                        risky.append((self.run, self.sub_run, count, gemroc))
                                        # print ("{} hits from Count= {},GEMROC ={}, FEB_couple={},Tcoarse[8:11]={}".format(self.tcoarse_dict[count][gemroc][FEB_group][page],count,gemroc,FEB_group,page))
                                    if (self.tcoarse_dict[count][gemroc][FEB_group][page] >= 32):
                                        sure_sat.append((self.run, self.sub_run, count, gemroc))

                                    if (self.tcoarse_dict[count][gemroc][FEB_group][page] > 32) :
                                        wierd.append((self.run, self.sub_run, count, gemroc))
                                        print (self.tcoarse_dict[count][gemroc][FEB_group][page] )
                                        print (f'Subrun: {self.sub_run},GEMROC {gemroc}, l1_count {count}')

                                    self.write_in_tree(count, gemroc, FEB_group, page)

                                except Exception as E:
                                    print ("\n--------------\nCount: {} wrong trigger time ({})".format(count,page))
                                    print (f'Trigger time {self.l1_trigger_dict[count]}')
                                    print (f'Trigger extremes {self.l1_trigger_extremes[count]}')
                                    print (f'Length key{self.length_dict[count].keys()}')
                                    print (E)

                                    pass

        self.generat_per_run_info(risky, sure_sat, wierd)
    def write_in_tree(self, count, gemroc, FEB_group, page):
        """
        Write one data in the tree
        :return:
        """
        # print (self.sub_run)
        self.tree_struct.run = int (self.run)
        self.tree_struct.subrun = int(self.sub_run)
        self.tree_struct.count = count
        self.tree_struct.gemroc= gemroc
        self.tree_struct.feb = FEB_group
        self.tree_struct.bucket = int(page,2)
        self.tree_struct.n_hit = self.tcoarse_dict[count][gemroc][FEB_group][page]

        self.to.Fill()

    def generat_per_run_info(self, risky, sure_sat, wierd):
        """
        Generate per count info and aggregate the data
        """
        count_list = []
        total_count_risky = 0
        total_count_sure = 0
        total_count_wierd=0

        for elem in risky:
            if elem[2] not in count_list:
                count_list.append(elem[2])
                total_count_risky += 1
            self.risky.append(elem)
        # print (count_list)
        count_list = []
        for elem in sure_sat:
            if elem[2] not in count_list:
                count_list.append(elem[2])
                total_count_sure += 1
            self.sure_sat.append(elem)

        count_list = []
        for elem in wierd:
            if elem[2] not in count_list:
                count_list.append(elem[2])
                total_count_wierd+= 1
            self.wierd.append(elem)

        self.total_count_risky += total_count_risky
        self.total_count_sure_sat += total_count_sure
        self.total_count_wierd += total_count_wierd

    def results(self):
        if self.total_trig!=0:
            print("------------------------------------------------")
            print(f'Number of triggers: {self.total_trig}')
            print(f'Number of risky FEBS: {len(self.risky)}')
            print(f'Number of sure_sat FEBS: {len(self.sure_sat)}')
            print(f'Number of wierd saturation (>32): {len(self.wierd)}')
            print(f'Number of trigger with risky situation: {self.total_count_risky}, {self.total_count_risky/self.total_trig*100:.2f}%')
            print(f'Number of trigger with sure_sat: {self.total_count_sure_sat}, {self.total_count_sure_sat/self.total_trig*100:.2f}%')



if __name__ == "__main__":
    num = sys.argv[1]
    runner = DPRAM_filling_check(num)
    runner.elab()
    #runner.elab_on_sub(80)
    runner.results()
