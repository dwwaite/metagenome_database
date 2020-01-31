# Cheat-sheet

A few quick commands to use while testing out the front-end of the database tool.

### Quick creation and value insertion

```bash
python metagenome_databaser.py create -d mock.hdb
python metagenome_databaser.py add -d mock.hdb --contig tests/contigs.fna \
    --prodigal_aa tests/genes_prod_aa.faa --prodigal_nt tests/genes_prod_nt.fna \
    --trna tests/genes_trna.fna --trna_contigs tests/contigs.fna \
    --rrna_ssu tests/genes_metaxa_ssu.fna \
    --coverage tests/coverage.txt --transcript tests/transcript.txt \
    --bin_file tests/bin_metadata.txt --bin_dir tests/bins
```

### View content of tables programmatically

```python
import sqlite3
db = sqlite3.connect('mock.hdb')

# Check tables are created
cursor = db.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(cursor.fetchall())

# Dump content of a table
cursor = db.cursor()
cursor.execute("SELECT * FROM contig;")
print(cursor.fetchall())
```