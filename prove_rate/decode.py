import binascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import ROOT
import sys
from ROOT import gROOT, AddressOf


class reader:
    def __init__(self, GEMROC_ID, MODE, RUN, SUBRUN):
        self.MODE = int(MODE)
        self.GEMROC_ID = int(GEMROC_ID)
        self.RUN = int(RUN)
        self.SUBRUN = int(SUBRUN)
        self.thr_scan_matrix = np.zeros((8, 64))  # Tiger,Channel
        self.thr_scan_frames = np.ones(8)
        self.thr_scan_rate = np.zeros((8, 64))
        self.warning_print_MAX=10
    def __del__(self):
        print("Done\n\n")

    def write_root(self, path):
        gROOT.ProcessLine('struct TreeStruct {\
				int runNo;\
				int layer;\
				int gemroc;\
				int tiger;\
				int channel;\
				int tac;\
				float tcoarse;\
				float tcoarse_10b;\
				float ecoarse;\
				float tfine;\
				float efine;\
				float charge_SH;\
				int count;\
				int count_ori;\
				int count_new;\
				int timestamp;\
				int l1ts_min_tcoarse;\
				int lasttigerframenum;\
				int local_l1_count;\
				int count_mismatch;\
				float delta_coarse;\
				int count_missing_trailer;\
				int subRunNo;\
                                int l1_framenum;\
                                int trailer_tiger;\
				};')
        from ROOT import TreeStruct

        rname = path.replace(".dat", ".root")

        # subRunNo = int(path.split("_")[7])
        subRunNo = self.SUBRUN

        # run = int(path.split("_")[6].split("/")[0])
        # print run

        rootFile = ROOT.TFile(rname, 'recreate')
        tree = ROOT.TTree('tree', '')
        mystruct = TreeStruct()

        for key in TreeStruct.__dict__.keys():
            if '__' not in key:
                formstring = '/F'
                if isinstance(mystruct.__getattribute__(key), int):
                    formstring = '/I'
                tree.Branch(key, AddressOf(mystruct, key), key + formstring)

        statinfo = os.stat(path)
        packet_header = -1
        packet_tailer = -1
        packet_udp = -1
        l1count = -1

        l1count_new = []
        lschannel = []
        lstac = []
        lstcoarse = []
        lstcoarse_10b = []
        lsecoarse = []
        lstfine = []
        lsefine = []
        lscharge_SH = []
        lstigerid = []
        lsl1ts_min_tcoarse = []
        lslasttigerframenum = []
        lscount_mismatch = []
        lsdelta_coarse = []
        l1timestamp = -1
        gemroc = -1
        l1framenum = -1
        trailer_tiger = -1

        pre_timestamp = 0
        pre_pretimestamp = 0

        tiger_framenum = -1
        prev_tiger_framenum = -1
        prev2_tiger_framenum = -1
        prev3_tiger_framenum = -1

        count_missing_trailer = 0

        hitcounter = 0
        max_hitcount = 1000000000
        flag_swap1 = False
        flag_swap2 = False
        firstPacket = True

        firstData = False
        print_debug = False

        with open(path, 'rb') as f:
            for i in range(0, statinfo.st_size // 8):
                data = f.read(8)
                # hexdata = binascii.hexlify(data)
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

                # for x in range(0, len(hexdata) - 1, 16):
                #     int_x = 0
                #     for b in range(7, 0, -1):
                #         hex_to_int = (int(str(hexdata[x + b * 2]), 16)) * 16 + int(str(hexdata[x + b * 2 + 1]), 16)
                #         int_x = (int_x + hex_to_int) << 8
                #
                #     hex_to_int = (int(str(hexdata[x]), 16)) * 16 + int(str(hexdata[x + 1]), 16)  # acr 2017-11-17 this should fix the problem
                #     int_x = (int_x + hex_to_int)

                ##############################################################################################
                ##																							##
                ##								TRIGGER-LESS DECODE											##
                ##																							##
                ##############################################################################################

                if (((int_x & 0xFF00000000000000) >> 59) == 0x00 and self.MODE == 0):
                    mystruct.runNo = self.RUN
                    mystruct.gemroc = self.GEMROC_ID
                    mystruct.tiger = (int_x >> 56) & 0x7
                    mystruct.channel = (int_x >> 48) & 0x3F
                    mystruct.tac = (int_x >> 46) & 0x3
                    mystruct.tcoarse = (int_x >> 30) & 0xFFFF
                    mystruct.ecoarse = (int_x >> 20) & 0x3FF
                    mystruct.tfine = (int_x >> 10) & 0x3FF
                    mystruct.efine = int_x & 0x3FF
                    mystruct.tcoarse_10b = (int_x >> 30) & 0x3FF

                    if (((int_x >> 20) & 0x3FF) - ((int_x >> 30) & 0x3FF)) > 0:
                        mystruct.delta_coarse = (((int_x >> 20) & 0x3FF) - ((int_x >> 30) & 0x3FF))
                    else:
                        mystruct.delta_coarse = (((int_x >> 20) & 0x3FF) - ((int_x >> 30) & 0x3FF)) + 1024
                    temp_ecoarse = mystruct.ecoarse
                    mystruct.charge_SH = int_x & 0x3FF

                    if (self.GEMROC_ID < 4):
                        mystruct.layer = 1
                    elif (self.GEMROC_ID > 3):
                        mystruct.layer = 2
                    if (self.GEMROC_ID > 11):
                        mystruct.layer = 0

                    tree.Fill()

                ##############################################################################################
                ##############################################################################################
                ##############################################################################################
                ##############################################################################################

                ##############################################################################################
                ##																							##
                ##								TRIGGER-MATCH DECODE										##
                ##																							##
                ##############################################################################################

                if (self.MODE == 1):
                    if (((int_x & 0xE000000000000000) >> 61) == 0x6):
                        # print "enter header"
                        packet_header = 1
                        LOCAL_L1_COUNT_31_6 = int_x >> 32 & 0x3FFFFFF
                        LOCAL_L1_COUNT_5_0 = int_x >> 24 & 0x3F
                        LOCAL_L1_COUNT = (LOCAL_L1_COUNT_31_6 << 6) + LOCAL_L1_COUNT_5_0
                        LOCAL_L1_TIMESTAMP = int_x & 0xFFFF
                        # print("local l1 count {}".format(LOCAL_L1_COUNT))
                        pre_pretimestamp = pre_timestamp
                        pre_timestamp = l1timestamp
                        # pre_l1count = l1count
                        l1count = LOCAL_L1_COUNT
                        l1timestamp = LOCAL_L1_TIMESTAMP

                        if firstData:
                            if print_debug:
                                print("WARNING: not able to record last tiger frame number of previous packet: no hits from TIGER 0-3 (L1_count={})!!!!!!!!!!!!!!!".format(LOCAL_L1_COUNT))

                        firstData = True  ## Header flags that next line will be first data word of the packet

                        if len(lschannel) > 0:
                            lschannel = []
                            lstac = []
                            lstcoarse = []
                            lsecoarse = []
                            lstfine = []
                            lsefine = []
                            lstcoarse_10b = []
                            lscharge_SH = []
                            lstigerid = []
                            lsl1ts_min_tcoarse = []
                            lslasttigerframenum = []
                            lscount_mismatch = []
                            l1count_new = []
                            lsdelta_coarse = []

                    if (((int_x & 0xC000000000000000) >> 62) == 0x0 and packet_header == 1 and packet_udp != 1):  ## DATA word
                        LOCAL_L1_TS_minus_TIGER_COARSE_TS = LOCAL_L1_TIMESTAMP - ((int_x >> 32) & 0xFFFF)
                        # print "enter DATA"
                        lstigerid.append((int_x >> 59) & 0x7)
                        lschannel.append((int_x >> 50) & 0x3F)
                        lstac.append((int_x >> 48) & 0x3)
                        lsecoarse.append((int_x >> 20) & 0x3FF)
                        lstfine.append((int_x >> 10) & 0x3FF)
                        lsefine.append(int_x & 0x3FF)
                        lslasttigerframenum.append((int_x >> 56) & 0x7)

                        lscharge_SH.append(int_x & 0x3FF)
                        temp_ecoarse = (int_x >> 20) & 0x3FF
                        lstcoarse_10b.append(((int_x >> 32) & 0x3FF))
                        temp_tcoarse = ((int_x >> 32) & 0x3FF)

                        tcoarse = (int_x >> 32) & 0xFFFF
                        ecoarse = (int_x >> 20) & 0x3FF
                        if (((int_x >> 20) & 0x3FF) - ((int_x >> 32) & 0x3FF)) > 0:
                            lsdelta_coarse.append(((int_x >> 20) & 0x3FF) - ((int_x >> 32) & 0x3FF))
                        else:
                            lsdelta_coarse.append(((int_x >> 20) & 0x3FF) - ((int_x >> 32) & 0x3FF) + 1024)

                        lstcoarse.append(tcoarse)

                        count_mismatch = 0

                        lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse
                        l1count_new_to_append = l1count

                        if int_x == 0:
                            print("WARNING: DATA with all zeros (subRun = {}, L1_count = {})".format(subRunNo, LOCAL_L1_COUNT))
                        else:
                            tiger_framenum = (int_x >> 56) & 0x7
                            if firstData:
                                prev3_tiger_framenum = prev2_tiger_framenum
                                prev2_tiger_framenum = prev_tiger_framenum
                                if (int_x >> 59) & 0x7 < 4:  ## store 2 previous tiger frame number (only from TIGER 0-3)
                                    prev_tiger_framenum = tiger_framenum
                                else:
                                    if print_debug:
                                        print("WARNING: not able to record last tiger frame number of this packet: no hits from TIGER 0-3 (L1_count={})!!!!!!!!!!!!!!!".format(LOCAL_L1_COUNT))

                                firstData = False

                            ########################################
                            ##         PACKETS MATCHING           ##
                            ########################################

                            ## Start from alignment of previous packet
                            if ((int_x >> 59) & 0x7 > 3):
                                if flag_swap1:
                                    temp_diff = pre_timestamp - tcoarse
                                elif flag_swap2:
                                    temp_diff = pre_pretimestamp - tcoarse
                                else:
                                    temp_diff = LOCAL_L1_TIMESTAMP - tcoarse
                            else:
                                temp_diff = lsl1ts_min_tcoarse_to_append  ## TIGER 0-3 always take the current l1ts

                            ## Find correct packet
                            ## performed only when lsl1ts_min_tcoarse is not inside the trigger window (roll-over taken into account)
                            if (not ((temp_diff > 1299 and temp_diff < 1567) or (temp_diff < -63960 and temp_diff > -64240))):

                                if firstPacket:  ## avoid packets correction for first packet
                                    pass  ## wrong entries in first packet should be discarded since they cannot be corrected
                                else:
                                    # print("Try SWAP 0")         					## try swap packets by 0
                                    temp_diff = LOCAL_L1_TIMESTAMP - tcoarse
                                    if ((temp_diff > 1299 and temp_diff < 1567) or (temp_diff < -63960 and temp_diff > -64240)):
                                        if flag_swap1 == True or flag_swap2 == True:
                                            print("SWAP 0 activated (L1_count={})".format(LOCAL_L1_COUNT))
                                        flag_swap1 = False
                                        flag_swap2 = False
                                    else:
                                        # print("Try SWAP 1")         					## try swap packets by 1
                                        temp_diff = pre_timestamp - tcoarse
                                        if ((temp_diff > 1299 and temp_diff < 1567) or (temp_diff < -63960 and temp_diff > -64240)):
                                            if flag_swap1 == False:
                                                print("SWAP 1 activated (L1_count={})".format(LOCAL_L1_COUNT))
                                            flag_swap1 = True
                                            flag_swap2 = False
                                        else:
                                            # print("Try SWAP 2")       					## try swap packets by 2
                                            temp_diff = pre_pretimestamp - tcoarse
                                            if ((temp_diff > 1299 and temp_diff < 1567) or (temp_diff < -63960 and temp_diff > -64240)):
                                                if flag_swap2 == False:
                                                    print("SWAP 2 activated (L1_count={})".format(LOCAL_L1_COUNT))
                                                flag_swap1 = False
                                                flag_swap2 = True
                                            else:
                                                if self.warning_print_MAX!=0:
                                                    print("WARNING: not able to correct packet (L1_count={}) !!!!!!!!!!!!!!!".format(LOCAL_L1_COUNT))
                                                    self.warning_print_MAX=self.warning_print_MAX-1

                            ## Apply packet correction to data of TIGER 4-7
                            if ((int_x >> 59) & 0x7 > 3):  ## correct packet for data of TIGER 4-7
                                if not (flag_swap1 or flag_swap2):  ## apply SWAP by 0 packet
                                    lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse  ## use l1ts of current packet
                                    l1count_new_to_append = l1count  ## use l1count of current packet
                                    count_mismatch = 0
                                elif flag_swap1:  ## apply SWAP by 1 packet
                                    lsl1ts_min_tcoarse_to_append = pre_timestamp - tcoarse  ## use l1ts of previous packet
                                    l1count_new_to_append = l1count - 1  ## use l1count of previous packet
                                    count_mismatch = 1
                                    if tiger_framenum != prev2_tiger_framenum:
                                        if print_debug:
                                            print("TIGER framecount not matched (SWAP1: L1_count={}: {} vs {}) !!!!!!!!!!!!!!!".format(LOCAL_L1_COUNT, tiger_framenum, prev2_tiger_framenum))
                                elif flag_swap2:  ## apply SWAP by 2 packets
                                    lsl1ts_min_tcoarse_to_append = pre_pretimestamp - tcoarse  ## use l1ts of 2 previous packet
                                    l1count_new_to_append = l1count - 2  ## use l1count of 2 previous packet
                                    count_mismatch = 2
                                    if tiger_framenum != prev3_tiger_framenum:
                                        if print_debug:
                                            print("TIGER framecount not matched (SWAP2: L1_count={}: {} vs {}) !!!!!!!!!!!!!!!".format(LOCAL_L1_COUNT, tiger_framenum, prev3_tiger_framenum))
                                else:
                                    print("Swap ERROR: a problem occurred in swap logic (subRun={}, L1_count={}) !!!!!!!!!!!!!!!".format(subRunNo, LOCAL_L1_COUNT))

                            ## Correct counters roll-over
                            if (lsl1ts_min_tcoarse_to_append < 0):
                                if ((int_x >> 59) & 0x7 > 3):
                                    if flag_swap1:
                                        lsl1ts_min_tcoarse_to_append = pre_timestamp - tcoarse + 2 ** 16
                                    if flag_swap2:
                                        lsl1ts_min_tcoarse_to_append = pre_pretimestamp - tcoarse + 2 ** 16
                                    if not (flag_swap1 or flag_swap2):
                                        lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse + 2 ** 16
                                else:
                                    lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse + 2 ** 16

                        #####################################################################################################################
                        #####################################################################################################################
                        #####################################################################################################################

                        lsl1ts_min_tcoarse.append(lsl1ts_min_tcoarse_to_append)
                        l1count_new.append(l1count_new_to_append)

                        lscount_mismatch.append(count_mismatch)

                    if (((int_x & 0xE000000000000000) >> 61) == 0x7):  ## TRAILER WORD --> sometimes is missing --> DO NOT USE
                        # print "enter trailer"
                        packet_tailer = 1
                        l1framenum = (int_x >> 37) & 0xFFFFFF
                        trailer_tiger = (int_x >> 27) & 0x7
                        gemroc = (int_x >> 32) & 0x1F

                    if (((int_x & 0xF000000000000000) >> 60) == 0x4):  ## UDP WORD --> used to flag end of packet
                        # print "enter UDP"
                        if packet_tailer == 0:
                            count_missing_trailer = count_missing_trailer + 1
                            print("WARNING: missing trailer word (subRun = {}, L1 count = {})!!!!!!!!!!!!!!!".format(subRunNo, LOCAL_L1_COUNT))
                        packet_udp = 1
                    # pre_udp_packet = udp_packet
                    # udp_packet = (((int_x >> 32)&0xFFFFF) + ((int_x >> 0) & 0xFFFFFFF))

                    if (packet_header == 1 and packet_udp == 1):  ## Fill ROOT file
                        for x in range(len(lstac)):
                            mystruct.channel = lschannel.pop()
                            mystruct.tac = lstac.pop()
                            mystruct.tcoarse = lstcoarse.pop()
                            mystruct.ecoarse = lsecoarse.pop()
                            mystruct.tfine = lstfine.pop()
                            mystruct.efine = lsefine.pop()
                            mystruct.tcoarse_10b = lstcoarse_10b.pop()
                            mystruct.charge_SH = lscharge_SH.pop()
                            mystruct.tiger = lstigerid.pop()
                            mystruct.l1ts_min_tcoarse = lsl1ts_min_tcoarse.pop()
                            mystruct.lasttigerframenum = lslasttigerframenum.pop()
                            mystruct.count_mismatch = lscount_mismatch.pop()
                            mystruct.count_new = l1count_new.pop()
                            mystruct.delta_coarse = lsdelta_coarse.pop()

                            mystruct.count_missing_trailer = count_missing_trailer
                            mystruct.local_l1_count = LOCAL_L1_COUNT
                            mystruct.count_ori = l1count
                            mystruct.count = mystruct.count_new
                            mystruct.timestamp = l1timestamp
                            mystruct.gemroc = self.GEMROC_ID  # previously was gemroc value from trailer word
                            mystruct.runNo = self.RUN
                            mystruct.subRunNo = subRunNo
                            mystruct.l1_framenum = l1framenum
                            mystruct.trailer_tiger = trailer_tiger

                            """
                            if(gemrocid<4):
                                mystruct.layer_id = 1
                            elif(gemrocid>3):
                                mystruct.layer_id = 2
                            if(gemrocid>11):
                                mystruct.layer_id = 0
                            """

                            if (self.GEMROC_ID < 4):
                                mystruct.layer = 1
                            elif (self.GEMROC_ID > 3):
                                mystruct.layer = 2
                            if (self.GEMROC_ID > 11):
                                mystruct.layer = 0

                            hitcounter = hitcounter + 1
                            if (hitcounter > max_hitcount):
                                sys.stderr.write("Increase max hit counter \n")
                                continue
                            # print "PRE FILL"
                            tree.Fill()
                        # print "POST FILL"
                        packet_header = 0
                        packet_tailer = 0
                        packet_udp = 0
                        firstPacket = False

        # print("Writing to ROOT file...")
        rootFile.Write()
        rootFile.Close()


# print("ROOT file written and closed.\n")


##############################################################################################
##																							##
##										MAIN												##
##																							##
##############################################################################################
if __name__ == "__main__":


    import sys

    print
    sys.argv

    if (len(sys.argv) == 4):
        filename = sys.argv[1]
        gemroc_id = sys.argv[2]
        mode = sys.argv[3]
        run = 1
    elif (len(sys.argv) == 5):
        filename = sys.argv[1]
        # run = int(filename.split("RUN_")[1].split("/")[0])
        subrun_number = int(os.path.basename(filename).split("_")[1])
        gemroc_id = sys.argv[2]
        mode = sys.argv[3]
        run = sys.argv[4]  # was wrong because Decode was given the wrong input from TER.sh
        print
        run
    # print subrun_number

    else:
        filename = input("Insert file name:")
        gemroc_id = input("Insert Gemroc id:")
        mode = input("MODE:(trigger_less:0 trigger_matched:1)")
        subrun_number=input("Subrun:")
        run = 1

    GEM5 = reader(gemroc_id, mode, run, subrun_number)
    print("Decoding: " + filename)
    GEM5.write_root(filename)
