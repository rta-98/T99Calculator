import os, sqlite3
import csv 
db_path = "/mnt/d/academic/tmp/test_1.db"
# To delete the file if it already exists 
# if os.path.exists(db_path):
#     os.remove(db_path)

conn = sqlite3.connect(db_path) # creates a new connection
cursor = conn.cursor() # creates a new cursor
cursor.execute("DROP TABLE IF EXISTS test_1") 
# --- End goal ---
# cursor.execute(""" 
# CREATE TABLE test_1 (
#     id INTEGER PRIMARY KEY, 
#     subst_iupac TEXT NOT NULL, 
#     subst_smiles TEXT UNIQUE,
#     temp_alpha_elim INTEGER NOT NULL,  
#     temp_beta_elim INTEGER NOT NULL, 
#     temp_alpha_clvg INTEGER NOT NULL,
#     init_prod_smarts TEXT, 
#     pdb_init_prod BLOB, 
#     head_id_smarts TEXT,
#     pdb_head_id BLOB
# )  
# """) 
cursor.execute(""" 
CREATE TABLE test_1 (
    id INTEGER PRIMARY KEY, 
    subst_iupac TEXT NOT NULL, 
    vtst_angstrom REAL, 
    temp_alpha_elim INTEGER,  
    temp_beta_elim INTEGER, 
    temp_alpha_clvg INTEGER,
    init_prod_smarts TEXT NOT NULL 
)  
""")

rows = [
    ("PFOA_linear", None, 700, 890, 940, "[#6](=[#8])(-[#9])-[*]"),
    ("PFOA_branched", None, 650, 810, 750, "[#6](=[#8])(-[#9])-[*]"),
    ("PFBA", None, 700, 890, 920, "[#6](=[#8])(-[#9])-[*]"), 
    ("HFPO-DA", None, 480, 860, 890, "[#6](=[#8])(-[*])-[#8]-[*]"),
    ("PFOS", 3.08, 640, 810, 610, "[#6](-[#6](-[#6](-[#6](-[#6](-[#6](-[#6](-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[#6](-[F])(-[F])-[F]"),
    ("PFBS", 3.04, 640, 810, 560, "[#6](-[#6](-[#6](-[#6](-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),
    ("FBSA", 3.08, 1030, 1060, 590, "[#6](-[#6](-[#6](-[#6](-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),  
    ("N-MeFBSAA", 3.18, None, None, 590, "[#6](-[#6](-[#6](-[#6](-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),
    ("2:3 FTCA", 2.79, 950, 1120, 880, "[#6](-[#6](-[#6](-[#6](-[#1])-[#1])(-[#1])-[#1])(-[F])-[F])(-[F])(-[F])-[F]"),
    ("3:2 FTCA", None, 920, 690, 950, "[#6](-[#6](-[#6](=[#6](-[#1])-[#1])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),
    ("2:2 FTS", None, 950, 710, 740, "[#6](-[#6](-[#6](=[#6](-[#1])-[#1])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),
    ("2:2 FTSA", None, 840, 650, 700, "[#6](-[#6](-[#6](=[#6](-[#1])-[#1])-[F])(-[F])-[F])(-[F])(-[F])-[F]"), 
    ("PFBOH", None, 470, 780, 840, "[#6](-[#6](-[#6](-[#6](-[F])=[#8])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[F]"), 
    ("2:2 FTOH", None, 1040, 770, 900, "[#6](-[#6](-[#6](=[#6](-[#1])-[#1])-[F])(-[F])-[F])(-[F])(-[F])-[F])"),
]

# Parametrized queries (?, ... ?) prevent SQL injection; <?> placeholders are used to bind data to the query.  
# INSERT opens a transaction which needs to be commited before changes are saved in the databsae. 
cursor.executemany("""
INSERT OR REPLACE INTO test_1
(subst_iupac, 
vtst_angstrom, 
temp_alpha_elim, 
temp_beta_elim, 
temp_alpha_clvg, 
init_prod_smarts) 
VALUES (?, ?, ?, ?, ?, ?) 
""", rows)

#|%%--%%| <eik9wM6aJk|NGyGXLqhrL>
conn.commit() # commits the previous transaction
for row in cursor.fetchall():
    print(row) 

#|%%--%%| <NGyGXLqhrL|YsYm80WTrg>
cursor.execute("SELECT * FROM test_1 ORDER BY id") 
rows = cursor.fetchall() 
with open("PFAS_thermal_data.csv", "w", newline="") as f: 
    write = csv.writer(f) 
    write.writerow([d[0] for d in cursor.description])
    write.writerows(rows)

#|%%--%%| <YsYm80WTrg|sLBpadKoVI>
conn.close() 
