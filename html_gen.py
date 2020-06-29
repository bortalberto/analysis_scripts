def wirte_html_run(run_number):
    with open("{}.html".format(run_number), 'w') as page:
        page.write("<!DOCTYPE html>\n \
            <html>\n \
            <head>\n \
            <title>RUN {0}</title>\n \
            </head>\n \
            <body>\n \
            <h1>Offline analysis result, RUN {0}</h1>\n ".format(run_number))

        page.write('<embed src="out/337.pdf" width="900" height="375"  type="application/pdf">')
        # page.write('<embed  src = "https://drive.google.com/viewerng/viewer?embedded = true & url = out/337.pdf " >')


        page.write("<p> Missing packets:</p>\n ")
        with open ("out/skipped_packets_run_{}".format(run_number), 'r') as misf:
            for line in misf.readlines():
                page.write("<p> {}</p> \n".format(line))

        page.write(" <img src='out/run_{}_sub_end.png'> \n".format(run_number))



        page.write("</body>\n \
            </html>\n ".format(run_number))
wirte_html_run(405)