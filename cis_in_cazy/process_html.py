

import os, subprocess, re

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
                    #sed processing on line
                    entry = line
                    entry = entry.rstrip()
                    entry = entry.lstrip()
                    entry = re.sub('<td>', '', entry)
                    entry = re.sub('</td>', '', entry)                    
                    entry = re.sub('<a href="http://www.cazy.org/', '', entry)
                    entry = re.sub('.html">', '', entry)
                    results[each_html].append(entry)
                    #results.append the next bit
                #process
            if "</table>" in line:
                print "done processing" + each_html
                break #done processing the table

for key, value in results.iteritems():
    for entry in value:
        print entry



            
for key, value in results.iteritems():
    line = key
    for each_field in value:
        line += each_field
    #print line
    #print len(value)
