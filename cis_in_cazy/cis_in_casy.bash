for all pdbs
curl cazy
curl "http://www.cazy.org/search?page=recherche&lang=en&recherche=1itx&tag=8" > 1itx.html


grep hits 1ubq.html
if "hits" appears:

parse these bits out to csv

<table class="search">
                <tr>
                        <th class="thsearch">Family</th>
                        <th class="thsearch">Kingdom</th>
                        <th class="thsearch">Organism</th>
                        <th class="thsearch">Protein Name</th>
                </tr>

                <tr>
                        <td><a href="http://www.cazy.org/GH18.html">
                                        GH18</a></td>
                        <td>Bacteria</td>
                        <td>Bacillus circulans WL-12</td>
                        <td>chitinase A1 (ChiA;ChiA1)</td>
                </tr>


        </table>



if hits does not appear:
report: not in cazy