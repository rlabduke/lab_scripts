

import os, subprocess

results = {}

for each_html in os.listdir("raw_html"):
    f = open("raw_html/"+each_html, 'r').read()
    if "hits" in f:
        hits = False #flag for if we hit the right part of the file yet
        for line in f:
            if "hits" in line:
                hits = True
                continue
            if hits:
                if "<td>" in line:
                    #results.append the next bit
                #process
            if "</table>" in line:
                break #done processing the table

#pseudocode for the rest of the project (do this in python)
# grep hits 1ubq.html
# if "hits" appears:

# parse these bits out to csv

# <table class="search">
#                 <tr>
#                         <th class="thsearch">Family</th>
#                         <th class="thsearch">Kingdom</th>
#                         <th class="thsearch">Organism</th>
#                         <th class="thsearch">Protein Name</th>
#                 </tr>

#                 <tr>
#                         <td><a href="http://www.cazy.org/GH18.html">
#                                         GH18</a></td>
#                         <td>Bacteria</td>
#                         <td>Bacillus circulans WL-12</td>
#                         <td>chitinase A1 (ChiA;ChiA1)</td>
#                 </tr>


#         </table>



# if hits does not appear:
# report: not in cazy
