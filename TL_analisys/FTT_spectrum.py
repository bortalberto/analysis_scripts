import ROOT as R
import array

def rescaleaxis(g,scale=1000/800):
    """This function rescales the x-axis on a TGraph."""
    N = g.GetN()
    x = g.GetX()
    for i in range(N):
        x[i] *= scale
    g.GetHistogram().Delete()
    g.SetHistogram(0)
    return

def extract_spectrum(tree,N_bins,ext_1,ext_2):
    h2 = R.TH1F("h2", 'rate Histogram', N_bins, ext_1, ext_2)
    for entryNum in range(0, tree.GetEntries()):
        tree.GetEntry(entryNum)
        timestamp = getattr(tree, "timestamp")
        if timestamp > ext_1 and timestamp < ext_2:
            h2.Fill(timestamp)
    h_trasf2 = h2.FFT(None, "MAG")
    gr2 = R.TGraph(h_trasf2)
    rescaleaxis(gr2 , 1 / (ext_2-ext_1))
    h_trasf_hz_2 = R.TH1F( "h_trasf_hz_2", 'rate histogram in Hz',N_bins , 0, 1/((ext_2-ext_1)/N_bins))
    x = array.array('d', [0])
    y = array.array('d', [0])
    for p in range(0, gr2 .GetN()):
        (gr2 .GetPoint(p, x, y))
        k = h_trasf_hz_2.FindFixBin(int(x[0]))
        h_trasf_hz_2.SetBinContent(k, int(y[0]))
    h_trasf_hz_2.SetXTitle("Rate [Hz]")
    # h_trasf_hz_2.GetXaxis().SetRangeUser(0, (ext_2-ext_1)/(N_bins*2))
    # spect = R.TSpectrum(10)
    # spect.Search(h_trasf_hz_2)
    # xs = spect.GetPositionX()
    # ys = spect.GetPositionY()
    # for i in range(0, 10):
    #     print(xs[i], ys[i])
    return h_trasf_hz_2