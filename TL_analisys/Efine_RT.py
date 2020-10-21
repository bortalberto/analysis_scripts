import ROOT as R
import array
import numpy as np
from math import log10
import pandas as pd
import datetime
import time
import deode_TL
import os
import FTT_spectrum as FFT_script
import socket
import tkinter as tk

gemroc=0
acq_time=1
running=False
c = R.TCanvas("c")
c.Draw()
c.SetLogy(1)
R.gStyle.SetOptStat(1)
def acquire_data(gemroc_num,run_time):
    HOST_IP="192.168.1.200"
    HOST_PORT=58880 + gemroc_num
    BUFSIZE=32000
    dataSock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    dataSock.settimeout(3)
    dataSock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    dataSock.bind((HOST_IP,HOST_PORT))
    data_list_tmp=[]
    start_time=time.time()
    while time.time()-start_time<run_time:
        data=dataSock.recv(BUFSIZE)
#         data="a"
        data_list_tmp.append(data)
    with open ("raw_data/efine_{}.raw".format(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), 'wb+') as savefile:
        for item in data_list_tmp:
            savefile.write(item)
            
            
decode_istance=deode_TL.reader("raw_data/")

def decode_data(filename,gemroc):
    out_name=("root_files/")+filename.split("raw_data/")[1].split(".")[0]+".root"
    decode_istance.write_root(filename,gemroc, 0, 0, out_name)
def get_last_file(mypath):
    onlyfiles = [f for f in os.listdir(mypath) if (os.path.isfile(os.path.join(mypath, f)) and 'missing_frames' not in f and "efine" in f)]
    onlyfiles.sort()
    return (os.path.join(mypath,onlyfiles[-1]))

def run_acq():
    while True:

        acquire_data(gemroc,acq_time)
        raw_f=get_last_file("raw_data/")
        decode_data(get_last_file("raw_data/"),gemroc)
        root_f=get_last_file("root_files/")

        f=R.TFile(get_last_file("root_files/"))
        f.tree.Draw("efine","")
        c.Update()
        os.remove(raw_f)
        os.remove(root_f)
        #root.update()

def status_run():
    running = True
    
def status_stop():
    running = False
    
#root = tk.Tk()
#frame = tk.Frame(root)
#frame.pack()

#button = tk.Button(frame, 
                   #text="Stop", 
                   #fg="red",
                   #command=status_stop)
#button.pack(side=tk.LEFT)

#slogan = tk.Button(frame,
                   #text="Start",
                   #command=status_run)
#slogan.pack(side=tk.LEFT)

run_acq()
