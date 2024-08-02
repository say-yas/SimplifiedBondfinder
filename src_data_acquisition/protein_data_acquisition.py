# 
# BSD 3-Clause License
# Copyright (c) 2023-2024 Sophia Bazzi and Sharareh Sayyad
# All rights reserved.

import Bio.PDB
from glob import glob
#import numpy as np
import pathlib
from pathlib import Path
import time
import os
import re
#import traceback
import math
from Bio.PDB import Selection, NeighborSearch, ShrakeRupley, calc_angle, calc_dihedral, MMCIFParser, PDBParser, parse_pdb_header
parser = PDBParser(QUIET=True)
mmcif_parser = MMCIFParser(QUIET=True)
from Bio.PDB.PDBExceptions import PDBConstructionException
import urllib.request
import PyPDF2
import warnings
warnings.filterwarnings("ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning)
from datetime import datetime
import platform
#import psutil

def process_atom(atom):
    """
    Extracts relevant information from the first target atom.

    Param:
        atom: The atom to be processed. 

    Returns:
        tuple: A tuple containing atom ID, residue name, residue number, chain ID, B-factor,
               occupancy, atom object, coordinates, and vector.
    """
    return (
        atom.get_id(),                             # Atom ID: eg. atom_S_info[0]
        atom.get_parent().get_resname(),           # Residue name
        atom.get_parent().get_id()[1],             # Residue number
        atom.get_parent().get_parent().get_id(),   # Chain ID 
        atom.get_bfactor(),                        # B-factor
        atom.get_occupancy(),                      # Occupancy
        atom,                                      # Atom object
        atom.coord,                                # Atom coordinates
        atom.get_vector()                          # Atom vector
    )

def process_neighbor_atom(neighbor):
    """
    Extracts relevant information from the atom adjacent to the first target atom.

    Param:
        neighbor: The neighboring atom to be processed.

    Returns:
        tuple: A tuple containing neighbor atom ID, residue name, residue number, occupancy,
               B-factor, atom ID, chain ID, vector, neighbor object, and coordinates.
    """
    return (
        neighbor.get_id(),                           # Neighbor atom ID
        neighbor.get_parent().get_resname(),         # Residue name
        neighbor.get_parent().get_id()[1],           # Residue number
        neighbor.get_occupancy(),                    # Occupancy
        neighbor.get_bfactor(),                      # B-factor
        neighbor.get_id(),                           # Redundant (Neighbor atom ID)
        neighbor.get_parent().get_parent().get_id(), # Chain ID
        neighbor.get_vector(),                       # Vector
        neighbor,                                    # Neighbor object
        neighbor.coord                               # Coordinates
    )
def process_structure(ii, structure, base_name):
    """
    Processes a protein structure and performs analysis.

    Param:
        ii (int): Index of the structure.
        structure: The structure to be processed.
        base_name (str): Base name of the structure file.
    """
    xray_model = structure[0]
    resolution = structure.header.get("resolution")
    method = structure.header.get("structure_method")

#    if resolution is not None and float(resolution) <= 2 and method == "X-RAY DIFFRACTION":
    if resolution is not None and float(resolution) <= 2:   
        print(resolution, method, structure)
        atoms = Bio.PDB.Selection.unfold_entities(xray_model, 'A')
        ns = Bio.PDB.NeighborSearch(atoms)

        for chain in xray_model.get_list():
            for residue in chain.get_list():
                for atom in residue.get_list():
                    if 'S' in atom.get_id() and atom.get_parent().get_resname() in Standard_residue:
                        atom_S_info = process_atom(atom)
                        close_atoms = [neighbor for neighbor in ns.search(atom.coord, 3.3) if neighbor != atom]

                        for neighbor in close_atoms:
                            if "N" in neighbor.get_id() and neighbor.get_id() not in blacklist and neighbor.get_parent().get_resname() in Standard_residue and neighbor.get_parent().get_id()[1] != atom_S_info[2]:
                                closeatom_N_info = process_neighbor_atom(neighbor)
                                process_CN_bond(atom_S_info, closeatom_N_info, chain, xray_model, ns, base_name, ii)

def process_CN_bond(atom_S_info, closeatom_N_info, chain,
                         xray_model, ns, base_name, ii):
    """
    Process the carbon atoms adjacent to the nitrogen atom to identify a CN bond
    """
    distance_SN = format(round(atom_S_info[6]-closeatom_N_info[8], 2), '.2f')
    print(distance_SN)
    CN_Neighbors = [neighbor for neighbor in ns.search(closeatom_N_info[9], 2.0) if "C" in neighbor.get_id() and atom_S_info[2] != closeatom_N_info[2] and neighbor.get_parent().get_id()[1] == closeatom_N_info[2] and neighbor.get_id() not in ["CD2", "CG"]]
    for CN_Neighbor in CN_Neighbors:
        Res_CN_num = CN_Neighbor.get_parent().get_id()[1]
        closeatom_C_info = process_neighbor_atom(CN_Neighbor)
        process_CS_bond(atom_S_info, closeatom_N_info, closeatom_C_info, chain,
                             xray_model, ns, base_name, ii, distance_SN)

def process_CS_bond(atom_S_info, closeatom_N_info, closeatom_C_info,
                         chain, xray_model, ns, base_name, ii, distance_SN):
    """
    Process the carbon atoms adjacent to the sulfur atom to identify a CS bond
    """
    # Search for sulfur (S) atom neighbors within a specified distance (2.50 Angstroms)
    XS_Neighbors = [neighbor for neighbor in ns.search(atom_S_info[7], 2.50) if neighbor.get_id() in S_C_ids and neighbor.get_parent().get_id()[1] == atom_S_info[2]]
    for CS_Neighbor in XS_Neighbors:
        Res_CS_num = CS_Neighbor.get_parent().get_id()[1]
        closeatom_CtoS_info = process_neighbor_atom(CS_Neighbor)
        structural_features_neighbour_residues(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
                             chain, xray_model, base_name, ii, distance_SN)

def structural_features_neighbour_residues(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
                         chain, xray_model, base_name, ii, distance_SN):
    
    """
    This function computes the distance, angle, and torsion angle between carbon, sulfur, and nitrogen atoms.
    Identifies residues within a 4-angstrom vicinity of the alpha carbon atoms of nitrogen- and sulfur-containing residues.
    """

    distance_CS = format(round(atom_S_info[6]-closeatom_CtoS_info[8], 2), '.2f')
    distance_CN = format(round(closeatom_N_info[8]- closeatom_C_info[8], 2), '.2f')
    angle_CSN = format(round(math.degrees(calc_angle(closeatom_CtoS_info[7], atom_S_info[8], closeatom_N_info[7])), 2), '.2f')
    angle_CNS = format(round(math.degrees(calc_angle(closeatom_C_info[7], closeatom_N_info[7], atom_S_info[8])), 2), '.2f')
    torsion_CS_NC =format(round(math.degrees(calc_dihedral(closeatom_CtoS_info[7], atom_S_info[8], closeatom_N_info[7], closeatom_C_info[7])), 2), '.2f')

    center_atomN = Selection.unfold_entities(chain[closeatom_N_info[2]], 'A')
    atom_list = [atom for atom in structure.get_atoms() if atom.name == 'CA']
    ns = NeighborSearch(atom_list)
    ngb_Res_toN = {res.get_resname() for center_atom in center_atomN for res in ns.search(center_atom.coord, 4, 'R')}
    
    center_atomS = Selection.unfold_entities(chain[atom_S_info[2]], 'A')
    ngb_Res_toS = {res.get_resname() for center_atom in center_atomS for res in ns.search(center_atom.coord, 4, 'R')}
    download_full_report(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
                          chain, xray_model, base_name, ii, distance_SN, distance_CS, distance_CN,
                            angle_CSN,angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS)

def download_full_report(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
                          chain, xray_model, base_name, ii, distance_SN, distance_CS, distance_CN,
                            angle_CSN, angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS):
    """
    This function downloads the full reports of each protein from the protein data bank.
    """
    pdb_id = base_name.split('_')[0]
    if not os.path.exists('Full_reports'):
       os.mkdir('Full_reports')
    output_folder = 'Full_reports'
    output_path = os.path.join(output_folder, f'{pdb_id}_full_validation.pdf')
    url = f'https://files.rcsb.org/pub/pdb/validation_reports/{pdb_id[1:3]}/{pdb_id}/{pdb_id}_full_validation.pdf'

    if pdb_id not in downloaded_files and atom_S_info[5] >= 0.8 and closeatom_N_info[3] >= 0.8:
       try:
          urllib.request.urlretrieve(url, output_path)
          print(f"Downloaded: {pdb_id}_full_validation.pdf")
          downloaded_files.append(pdb_id)
       except Exception as exc:
           print(f"Failed to download {pdb_id}_full_validation.pdf. Error: {exc}")
    if atom_S_info[5] >= 0.8 and closeatom_N_info[3] >= 0.8:
       RSRZ(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
             chain, xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
               angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, output_path)
    else:
      print('atom_occupation_below_threshold')
    

def RSRZ(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info,
          chain, xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
            angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, output_path):
    """
    This function reads Full_reports file in pdf format containing the Real Space R-factor Z-score (RSRZ) data and identifes residues that are RSRZ outliers
    Records outlier residues in an output file.
    If no outliers are found, calls the `sasa` function for further processing.
    """
    Scompletelist = set()
    Ncompletelist = set()
    found_match = False
    print(f'processing {closeatom_N_info[1]}_{closeatom_N_info[0]}_{closeatom_N_info[2]}_{closeatom_N_info[6]}, {atom_S_info[1]}_{atom_S_info[0]}_{atom_S_info[2]}_{atom_S_info[3]} for pdb_id: {base_name.replace("_final", "")}')
    with open(output_path, 'rb') as outlier_file:
         pdf_reader = PyPDF2.PdfFileReader(outlier_file, strict=False)
         pages = pdf_reader.numPages
         for jj in range(pages):
             pageObj = pdf_reader.getPage(jj)
             text = pageObj.extractText().split("\n")
             if all(keyword in text for keyword in ["Mol" and "Chain" and "Res" and "Type" and "RSRZ"]):
                for idx, item in enumerate(text):
                    if item == str(atom_S_info[2]):
                       before = text[idx - 1] if idx - 1 >= 0 else ""
                       after = text[idx + 1] if idx + 1 < len(text) else ""
                       resSnum, chainS, resS  = item, before, after
                       Scomplete = f"{resS}_{resSnum}_{chainS}"
                       Scompletelist.add(Scomplete)
                    elif item == str(closeatom_N_info[2]):
                       before = text[idx - 1] if idx - 1 >= 0 else ""
                       after = text[idx + 1] if idx + 1 < len(text) else ""
                       resNnum, chainN, resN  = item, before, after
                       Ncomplete = f"{resN}_{resNnum}_{chainN}"
                       Ncompletelist.add(Ncomplete)
                
                # If any data has been found, set flag to true
                if Scompletelist or Ncompletelist:
                    found_match = True
    
    outliers = Scompletelist | Ncompletelist

    existing_lines = set()
    if os.path.exists('RSRZ_outliers.txt'):
        with open('RSRZ_outliers.txt', 'r') as existing_file:
            existing_lines = set(existing_file.readlines())
    
    with open('RSRZ_outliers.txt', 'a') as outliers_file:
        for outlier in outliers:
            line = f"the residue {outlier} in {base_name.replace('_final', '')} is an RSRZ outlier\n"
            if line not in existing_lines:
               outliers_file.write(f"the residue {outlier} in {base_name.replace('_final', '')} is an RSRZ outlier\n")
    # If no match was found across all pages, call the sasa function
    if not found_match:
        sasa(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain, xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN, angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, output_path)

def sasa(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
          xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
            angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, output_path):
    
    """
    This function computes the solvent accessible surface area (SASA) of the sulfur and nitrogen atoms.
    """

    SASA_Res_S = None
    SASA_Res_N = None
    try:
      sr = ShrakeRupley()
      sr.compute(xray_model, level="R")
      SASA_Res_N = round(xray_model[closeatom_N_info[6]][closeatom_N_info[2]].sasa, 2)
      SASA_Res_S = round(xray_model[atom_S_info[3]][atom_S_info[2]].sasa, 2)

      SASA_atom_N = round(xray_model[closeatom_N_info[6]][closeatom_N_info[2]][closeatom_N_info[0]].sasa, 2)
      SASA_atom_S = round(xray_model[atom_S_info[3]][atom_S_info[2]][atom_S_info[0]].sasa, 2)
      
      download_BDB(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                    xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                      angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                        SASA_Res_S, SASA_atom_S, output_path)
    except Exception as exc:
        print(exc)
        pass

def download_BDB(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                  xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                    angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                      SASA_Res_S, SASA_atom_S, output_path):
    
    """
    This function downloads files from BDB.
    """

    pdb_id = base_name.split('_')[0]
    if not os.path.exists('BDB'):
       os.mkdir('BDB')
    output_folderBDB = 'BDB'
    output_path_BDB = os.path.join(output_folderBDB, f'{pdb_id}.bdb')
    url = f'https://www3.cmbi.umcn.nl/bdb/download/{pdb_id}/'
    
    if pdb_id not in downloaded_BDB_files:
       try:
          urllib.request.urlretrieve(url, output_path_BDB)
          print(f"Downloaded: {pdb_id}.bdb")
          downloaded_BDB_files.append(pdb_id)
       except Exception as exc:
          print(f"Failed to download {pdb_id}.bdb. Error: {exc}")
          pass
    process_BDB(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                 xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                   angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                     SASA_Res_S, SASA_atom_S, output_path_BDB, output_path)

def process_BDB(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                 xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                   angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                     SASA_Res_S, SASA_atom_S, output_path_BDB, output_path):
    
    """
    This function processes the BDB to extract the B-factors of the sulfur and nitrogen atoms.
    """

    bfac_S_BDB = None
    bfac_N_BDB = None
    try:
        parser = PDBParser(QUIET=True)
        structureBB = parser.get_structure(ii, output_path_BDB)
        xray_model = structureBB[0]
        atoms = Bio.PDB.Selection.unfold_entities(xray_model, 'A')
        ns = Bio.PDB.NeighborSearch(atoms)
        print('BDB', structureBB)

        for chain in xray_model.get_list():
            for residue in chain.get_list():
                for atom in residue.get_list():
                    if 'S' in atom.get_id() and atom.get_parent().get_resname() in Standard_residue and atom.get_parent().get_parent().get_id() == atom_S_info[3]:
                        bfac_S_BDB = atom.get_bfactor() 
                        close_atoms = [neighbor for neighbor in ns.search(atom.coord, 3.3) if neighbor != atom]

                        for neighbor in close_atoms:
                            if "N" in neighbor.get_id() and neighbor.get_id() not in blacklist and neighbor.get_parent().get_resname() in Standard_residue and neighbor.get_parent().get_id()[1] != atom_S_info[2]:
                                closeatom_N = neighbor
                                bfac_N_BDB = closeatom_N.get_bfactor()

        too_close_contacts(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                            xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                              angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                                SASA_Res_S, SASA_atom_S, output_path_BDB, output_path, bfac_S_BDB, bfac_N_BDB)

    except Exception as exc:
        print("Exception:", exc)
        pass 
        too_close_contacts(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                            xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                              angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                                SASA_Res_S, SASA_atom_S, output_path_BDB, output_path, bfac_S_BDB, bfac_N_BDB)

def too_close_contacts(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
                        xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
                          angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                            SASA_Res_S, SASA_atom_S, output_path_BDB, output_path, bfac_S_BDB, bfac_N_BDB):
    
    """
    This fuction collects the close contacts between the sulfur and nitrogen atoms from Full_reports file in pdf format.
    """
    
    search_term1 = f'{atom_S_info[3]}:{atom_S_info[2]}:{atom_S_info[1]}:{atom_S_info[0]}'
    search_term2 = f'{closeatom_N_info[6]}:{closeatom_N_info[2]}:{closeatom_N_info[1]}:{closeatom_N_info[0]}'
    close_contacts = None
    with open(output_path, 'rb') as outlier_file:
         pdf_reader = PyPDF2.PdfFileReader(outlier_file)
         header_pattern = re.compile(r'\b\s*Atom-1\s*\s*Atom-2\s*\s*Interatomic\s*\b', re.IGNORECASE)
         found_term1 = False
         found_term2 = False
         for page_num in range(len(pdf_reader.pages)):
             page = pdf_reader.getPage(page_num)
             page_text = page.extractText()
             match = header_pattern.search(page_text)

             if match:
                table_page = pdf_reader.getPage(page_num)
                table_text = table_page.extractText()
                rows = [row.split() for row in table_text.split('\n') if row.strip()]
                for row in rows:
                    if any(search_term1 in element for element in row):
                       found_term1 = True
                    if any(search_term2 in element for element in row):
                       found_term2 = True
                    if found_term1 and found_term2:
                       close_contacts = True
                       break
                       return close_contacts
                       print(True)
   
    writing(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
             xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
               angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                 SASA_Res_S, SASA_atom_S, output_path_BDB, output_path, bfac_S_BDB,
                   bfac_N_BDB, found_term1, found_term2, close_contacts)

def writing(atom_S_info, closeatom_N_info, closeatom_C_info, closeatom_CtoS_info, chain,
             xray_model, base_name, ii, distance_SN, distance_CS, distance_CN, angle_CSN,
               angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, SASA_Res_N, SASA_atom_N,
                 SASA_Res_S, SASA_atom_S, output_path_BDB, output_path, 
                 bfac_S_BDB, bfac_N_BDB, found_term1, found_term2, close_contacts):

    if close_contacts:
       close_contacts = "Close"
    else:
       close_contacts = "Not-close"

    fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\n'.format(
    atom_S_info[3], atom_S_info[1], atom_S_info[0], atom_S_info[2], atom_S_info[4], bfac_S_BDB, atom_S_info[5], SASA_Res_S, SASA_atom_S,
    closeatom_N_info[6], closeatom_N_info[1], closeatom_N_info[0], closeatom_N_info[2], closeatom_N_info[4], bfac_N_BDB,
    closeatom_N_info[3], SASA_Res_N, SASA_atom_N, closeatom_C_info[0], closeatom_C_info[4], closeatom_C_info[3],
    closeatom_C_info[1], closeatom_C_info[0], closeatom_C_info[4], closeatom_C_info[3], distance_SN,
    distance_CN, distance_CS, angle_CSN, angle_CNS, torsion_CS_NC, ngb_Res_toN, ngb_Res_toS, close_contacts, f"{base_name.replace('_final', '')}_{atom_S_info[3]}_{atom_S_info[1]}_{atom_S_info[0]}_{atom_S_info[2]}_{closeatom_N_info[6]}_{closeatom_N_info[1]}_{closeatom_N_info[0]}_{closeatom_N_info[2]}.pdb"))

# The code starts here
BaseDir = os.getcwd()
def format_time(seconds):
    hours, rem = divmod(seconds, 3600)
    minutes, rem = divmod(rem, 60)
    seconds, milliseconds = divmod(rem, 1)
    milliseconds = int(milliseconds * 1000)
    return f"{int(hours)}h {int(minutes)}m {int(seconds)}s {milliseconds}ms"

D = Path('.')
BDB_path = Path('./BDB/')
blacklist = ['ZN', 'H4N']
S_C_ids = ['CB', 'CG']
Standard_residue = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'PYL', 'SEC']
#
with open('standard_restest.txt', 'a', 1) as fout:
    fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\n'.format(
        "Chain_S", "Res_S", "S_id", "Res_S_num", "bfac_S", "bfac_S_BDB", "occ_S", "SASA_S", "SASA_atom_S", "Chain_N", "Res_N", "N_id", "Res_N_num",
        "bfac_N", "bfac_N_BDB", "occ_N", "SASA_N", "SASA_atom_N", "CN_id", "bfac_CN", "occu_CN", "Res_CS", "CS_id", "bfac_CS", "occu_CS",
        "distance_SN", "distanceCN", "distanceCS", "angle_CSN", "angle_CNS", "torsion_CS_NC", "ngb_Res_toN",
        "ngb_Res_toS", "Close_contact", "filename"))

    total_time = 0
    file_count = 0
    downloaded_files = []
    downloaded_BDB_files = []
    for root, dirs, files in os.walk(D):
        for ii, filename in enumerate(files, 1):
            start_time = time.time()
            if filename.endswith("_final.cif"):
               base_name = Path(filename).stem
               file_path = os.path.join(root, filename)
               print(file_path, base_name)

               # Increment the file count for each matching file
               file_count += 1

               try:
                 parser = MMCIFParser()
                 structure = parser.get_structure(ii, file_path)
                 process_structure(ii, structure, base_name)
                # print('parser', structure)
               except Exception as exc:
                  print("Exception:", exc)
                  traceback.print_exc()
                  pass
     
               finally:
                   end_time = time.time()
                   execution_time = end_time - start_time
                   total_time += execution_time
                   average_time = total_time / file_count
                   print(f"Processing Time for {filename}: {format_time(execution_time)}")
current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print("Date:", current_date)
print(f"Total Time: for {file_count} files", format_time(total_time))
print(f"Average Time per iteration: {format_time(average_time)}")
