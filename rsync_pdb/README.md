The PDB rsyncing script is maintained at ftp://ftp.rcsb.org/pub/pdb/software/rsyncPDB.sh.  Open it up, fix some paths, and run it to do the rsync.

A list of all the files can be generated with:

> find /home/smlewis/whole_PDB_as_PDB . -name "*ent.gz"


/home/smlewis/whole_PDB_as_PDB/w0/pdb3w0q.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb1w0g.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb2w0a.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb1w0s.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb3w0g.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb1w02.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb3w0j.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb3w08.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb3w0o.ent.gz
/home/smlewis/whole_PDB_as_PDB/w0/pdb1w04.ent.gz
...


Pipe to a file to save it, or pipe through "sort -R" if you want them in random order for some reason.

Note the odd file format: PDB supplies 1UBQ as "pdb1ubq.ent.gz", not 1UBQ.pdb.gz as might have been expected.  Also note the directory structure: the PDBs are organized into directores by the MIDDLE TWO characters of their 4-character label.