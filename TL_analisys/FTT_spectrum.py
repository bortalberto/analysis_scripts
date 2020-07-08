import ROOT as R
import array
from math import log10
def rescaleaxis(g,scale=1000/800):
    """This function rescales the x-axis on a TGraph."""
    N = g.GetN()
    x = g.GetX()
    for i in range(N):
        x[i] *= scale
    g.GetHistogram().Delete()
    g.SetHistogram(0)
    return

# def rescaleaxis2(g,scale=1000/800):
#     """This function rescales the x-axis on a TGraph."""


def extract_spectrum(tree,rate_max,ext_1,ext_2,remove_peaks=False):
    N_bins=int((ext_2-ext_1)*rate_max)
    rate_real=N_bins/(ext_2-ext_1)
    h2 = R.TH1F("h2", 'rate Histogram', N_bins, ext_1, ext_2)
    tree.Draw("timestamp>>h2")
    if remove_peaks:
        for num in range (0, N_bins):
            if h2.GetBinContent(num)>300:
                h2.SetBinContent(num,h2.GetBinContent(num-1))
    h_trasf2 = h2.FFT(None, "MAG")
    gr2 = R.TGraph(h_trasf2)
    rescaleaxis(gr2 , rate_real/N_bins)
    h_trasf_hz_2 = R.TH1F( "h_trasf_hz_2", 'rate histogram in Hz',int(N_bins/2 ), 0, rate_real/2)
    x = array.array('d', [0])
    y = array.array('d', [0])
    for p in range(0, gr2 .GetN()):
        (gr2 .GetPoint(p, x, y))
        if x[0]<rate_real/2:
            k = h_trasf_hz_2.FindFixBin((x[0]))
            # print (x,y, k)
            h_trasf_hz_2.SetBinContent(k, (y[0]))
    h_trasf_hz_2.SetXTitle("Rate [Hz]")
    # h_trasf_hz_2.GetXaxis().SetRangeUser(0, (ext_2-ext_1)/(N_bins*2))
    # spect = R.TSpectrum(10)
    # spect.Search(h_trasf_hz_2)
    # xs = spect.GetPositionX()
    # ys = spect.GetPositionY()
    # for i in range(0, 10):
    #     print(xs[i], ys[i])
    h_spec = R.TH1F( "h_spec", 'PSD vs Hz',int(N_bins/2) , 0, rate_real/2)
    for n in range(0, h_trasf_hz_2.GetNbinsX()):
        h_spec.SetBinContent(n, (h_trasf_hz_2.GetBinContent(n) ** 2 / N_bins))
    max=h_spec.GetMaximum()
    for n in range(0, h_spec.GetNbinsX()):
        if (h_spec.GetBinContent(n) / max) !=0:
            h_spec.SetBinContent(n, 10*log10((h_spec.GetBinContent(n) / max)))
        else:
            h_spec.SetBinContent(n,-100)
    return (h_trasf_hz_2.Clone(),h_spec.Clone(),h_trasf2.Clone())



#
# def extract_spectrum_2(tree,rate_max,ext_1,ext_2):
#     N_bins=int((ext_2-ext_1)*rate_max)
#     h2 = R.TH1F("h2", 'rate Histogram', N_bins, ext_1, ext_2)
#     tree.Draw("timestamp>>h2")
#     h_trasf2 = h2.FFT(None, "MAG")
#
#
#     h_trasf_hz_2.SetXTitle("Rate [Hz]")
#     # h_trasf_hz_2.GetXaxis().SetRangeUser(0, (ext_2-ext_1)/(N_bins*2))
#     # spect = R.TSpectrum(10)
#     # spect.Search(h_trasf_hz_2)
#     # xs = spect.GetPositionX()
#     # ys = spect.GetPositionY()
#     # for i in range(0, 10):
#     #     print(xs[i], ys[i])
#     h_spec = R.TH1F( "h_spec", 'PSD vs Hz',N_bins , 0, 1/((ext_2-ext_1)/N_bins))
#     for n in range(0, h_trasf_hz_2.GetNbinsX()):
#         h_spec.SetBinContent(n, (h_trasf_hz_2.GetBinContent(n) ** 2 / N_bins))
#     max=h_spec.GetMaximum()
#     # for n in range(0, h_spec.GetNbinsX()):
#     #     if (h_trasf_hz_2.GetBinContent(n) / max) !=0:
#     #         h_spec.SetBinContent(n, 10*log10((h_trasf_hz_2.GetBinContent(n) / max)))
#     #     else:
#     #         h_spec.SetBinContent(n,0)
#     return (h_trasf_hz_2,h_spec)