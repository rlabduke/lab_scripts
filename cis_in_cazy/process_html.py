

import os, subprocess

results = {}

for each_html in os.listdir("raw_html"):
    f = open("raw_html/"+each_html, 'r').read()
    #flines = f.splitlines()
    if "hits" in f:
        hits = False #flag for if we hit the right part of the file yet
        #print each_html
        results[each_html] = []
        for line in f.splitlines():
            #print line
            if "hits" in line:
                hits = True
                #print "hitsflag"
                continue
            if hits:
                if "<td>" in line:
                    #print line
                    results[each_html].append(line)
                    #results.append the next bit
                #process
            if "</table>" in line:
                print "done processing" + each_html
                break #done processing the table

print results
