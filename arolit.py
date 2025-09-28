import subprocess
import pymupdf4llm
import sqlite3
import pandas as pd
import os
import re
import time
import csv
from timeoutwrapper import timeout

#%% Main 
class arolit:

    # Initialize class
    def __init__(self):
        self.doi_dict = {}
        self.PMID_dict = {}
        self.primers_dict = []
        pass

    def fetch_pubmed_medline(self, query: str, output_file: str):
        """
        Fetch MEDLINE records from PubMed using Entrez Direct (via WSL).

        Simulates:
            wsl bash -c "esearch -db pubmed -query 'query' | efetch -format medline"

        Args:
            query (str): PubMed search query.
            output_file (str): Path to save the MEDLINE result.
        """
        try:
            # Build full command as a string for bash
            full_cmd = f"esearch -db pubmed -query \"{query}\" | efetch -format medline"

            # Run it inside WSL bash
            result = subprocess.run(
                ["wsl", "bash", "-c", full_cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            if result.returncode != 0:
                print("Command failed:", result.stderr.decode())
                return

            # Write output to file
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(result.stdout.decode())

            return(f"Successfully fetched articles.")

        except Exception as e:
            return(f"Error running command in WSL: {e}")
        
            
    def fetch_doi_dict(self, nbib_file):
        self.doi_dict = {}
        self.PMID_dict = {} # some articles don't have DOIs so we can try to save PMID instead (always has PMID bc it was fetched from PubMed db)
        pmid = None # so it exists outside of the loop for safety
        has_doi = False # initialize condition
        
        with open(nbib_file, 'r', encoding='utf-8') as file:
            for line in file:
                if line.startswith('PMID- '): 
                    pmid = line.strip().split('- ')[1]
                    has_doi = False # reset because new article
                    
                elif line.startswith('AID - ') and '[doi]' in line:
                    match = re.search(r'10\.\d{4,9}/\S+', line)
                    if match:
                        doi = match.group(0)
                        self.doi_dict.setdefault(f"https://doi.org/{doi}", [])
                        has_doi = True # mark article as having DOI
                        
                elif line.strip() == '': # spots the empty line inbetween entries
                    if pmid and not has_doi:
                        self.PMID_dict.setdefault(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}", [])
                    pmid = None
                    has_doi = False # reset back to default second time for safe measures


        return self.doi_dict, self.PMID_dict

    @timeout(300) # 5 minute timeout
    def extract_text_from_pdf(self, file_path): 
        text = pymupdf4llm.to_markdown(file_path)
        return text
    
    def oligo_search_new(self,article):
        pattern = re.compile(
            r'''(?ix)
            (?<![A-Za-z0-9])                                  # left boundary (no letter/number just before)
            (?:                                                # start sequence
            [ACGTURKSMYWBHNDV]                               # 1st base
            (?:[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]*  # separators: space/tab/newline or PDF dash variants
            [ACGTURKSMYWBHNDV]                            # next base
            ){16,}                                           # => total ≥17 bases
            )
            (?![A-Za-z0-9])                                    # right boundary (no letter/number just after)
        ''')
        clean = lambda s: re.sub(r'[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]', '', s).upper()

        hits = [clean(m.group(0)) for m in pattern.finditer(article)]
        # keep only unique, length-validated primers (safety net)
        seen = set()
        hits = [h for h in hits if len(h) >= 17 and not (h in seen or seen.add(h))]
        return hits

    def oligo_search_new2(self, article):
        pattern = re.compile(
            r'''(?ix)
            (?<![A-Za-z0-9])                                  # left boundary (no letter/number just before)
            (?:                                                # start sequence
                [ACGTURKSMYWBHNDV/()]                               # 1st base
                (?:
                    [\s\-\u2010\u2011\u2012\u2013\u2014\u2212]*  # separators: space/tab/newline or PDF dash variants
                    [ACGTURKSMYWBHNDV/()]*                            # next base
                )                                           # => total ≥17 bases
            )
            (?![A-Za-z0-9])                                    # right boundary (no letter/number just after)
        ''')

        clean = lambda s: re.sub(r'[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]', '', s).upper()
        hits = [clean(m.group(0)) for m in pattern.finditer(article)]
        hits = [hit.replace(' ','') for hit in hits]
        hits = [hit.strip('()') for hit in hits]
        # keep only unique, length-validated primers (safety net)
        seen = set()
        hits = [h for h in hits if len(h) >= 17 and not (h in seen or seen.add(h))]
        return hits

    def oligos_to_csv_mult(self, files, names, output_csv):

        # files = [f for f in os.listdir(folder_path) if f.endswith('.pdf')]
        data = []

        for file, name in zip(files, names):
            # path = os.path.join(folder_path, file)
            print(f'Processing {name}...')
            
            doi = os.path.splitext(name)[0] # Placeholder, objetivo é que o titulo seja o DOI. Esta função divide o nome numa tupla com antes e depois da extensão daí o [0]
            
            try:
                # Extrair texto do PDF
                text = self.extract_text_from_pdf(file)

                oligos = self.oligo_search_new(text)
                
                oligos_2 = self.oligo_search_new2(text)

                for sequence in oligos:
                    data.append([sequence, doi])

                seen_sequences = {row[0] for row in data}

                for sequence in oligos_2:
                    if sequence not in seen_sequences:
                        data.append([sequence, doi])
                        seen_sequences.add(sequence)

                self.primers_dict = data

            except Exception as e:
                return(f'Skipping {file} due to error: {e}')
        
        with open(output_csv, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(['Sequence', 'Source'])
            writer.writerows(data)
        
        return 'Primers successfully parsed to CSV.'

    def csv_to_sql(self, csv_path, database):
        conn = sqlite3.connect(database)
        cursor = conn.cursor()
        
        # Create tables if they don't exist
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sequences (
                ID INTEGER PRIMARY KEY AUTOINCREMENT,
                Sequence TEXT UNIQUE
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS [references] (
                ID INTEGER PRIMARY KEY AUTOINCREMENT,
                Reference TEXT UNIQUE
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sequence_references (
                sequence_id INTEGER,
                reference_id INTEGER,
                FOREIGN KEY (sequence_id) REFERENCES sequences(ID),
                FOREIGN KEY (reference_id) REFERENCES [references](ID),
                PRIMARY KEY (sequence_id, reference_id)
            )
        """)
        # Reset das tabelas para evitar problemas com o autoincrement
        cursor.execute("DELETE FROM sequence_references")  
        cursor.execute("DELETE FROM sequences")            
        cursor.execute("DELETE FROM \"references\"")
        
        cursor.execute("DELETE FROM sqlite_sequence WHERE name='sequences'")
        cursor.execute("DELETE FROM sqlite_sequence WHERE name='references'")

        # Carrega o ficheiro CSV em forma de dataframe, onde cada coluna do CSV é uma coluna da dataframe
        df = pd.read_csv(csv_path)

        # Insere sequências na base de dados
        for sequence in df['Sequence'].unique(): # Extrai todas as sequências únicas da coluna Sequences
            cursor.execute("INSERT OR IGNORE INTO sequences (sequence) VALUES (?)", (sequence,)) # Ignora o Insert se a sequência já existir

        # Insere referências na base de dados
        for reference in df['Source'].unique(): # Igual ao de cima
            cursor.execute("INSERT OR IGNORE INTO \"references\" (reference) VALUES (?)", (reference,)) # Usar \"references\" pois é uma keyword reservada no SQL

        # Junta sequências com as respetivas referências
        for _, row in df.iterrows(): # Itera sobre cada linha. Como dá return a um par Index-Conteúdo, ignoramos o index com '_'
            sequence = row['Sequence']
            reference = row['Source']
            cursor.execute("""
                INSERT OR IGNORE INTO sequence_references (sequence_id, reference_id)
                VALUES (
                    (SELECT id FROM sequences WHERE sequence = ?),
                    (SELECT id FROM "references" WHERE reference = ?)
                )
            """, (sequence, reference))

        # Commit changes and close connection
        conn.commit()
        conn.close()
        return "Database created successfully."

    def export_for_score(self, input_db, output_csv):
        # Creates CSV file with indexed sequences and the articles they show up in

        # Connect to the database
        ''' Database used for input should have 3 tables: Sequences (contains an
        ID column and the corresponding sequence), References (contains an ID
        column and the corresponding reference), and Sequences_references 
        (contains a column with the sequence ID from Sequences and the reference
        ID from References) '''
        conn = sqlite3.connect(input_db)
        cursor = conn.cursor()
            
        # Run the query
        query = """
        SELECT s.ID, s.Sequence, GROUP_CONCAT(r.Reference, ' | ') as Ref_list
        FROM sequences s
        LEFT JOIN sequence_references sr ON s.ID = sr.sequence_id
        LEFT JOIN "references" r ON sr.reference_id = r.ID
        GROUP BY s.ID, s.Sequence;
        """
        data = cursor.execute(query).fetchall()
            
        # Close connection
        conn.close()
            
        # Convert to DataFrame and format Sequence ID
        df = pd.DataFrame(data, columns=["Sequence ID", "Sequence", "References"])
        df["Sequence ID"] = df["Sequence ID"].apply(lambda x: f"Ebola{x}")
            
        # Convert DOIs into full links
        def format_doi(ref_string):
            if ref_string is None:
                return ""
            doi_list = ref_string.split(" | ")
            # Change this as needed depending on how your DOI was saved as
            formatted_dois = [f"https://doi.org/{doi[0:7]}/{doi[7:]}" for doi in doi_list]
            return " | ".join(formatted_dois)
            
        # Apply formatting to the References column
        df["References"] = df["References"].apply(format_doi)
            
        # Calculate primer length
        df.insert(2, "Primer length", [len(seq) for seq in df["Sequence"].tolist()])
            
        # Save to CSV
        df.to_csv(output_csv, index=False)
            
        # Save to primers_dict
        self.primers_dict = df.to_dict(orient='records')
            
        return(f"CSV file saved to {output_csv}")