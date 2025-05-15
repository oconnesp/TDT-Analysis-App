import re
import numpy as np
from datetime import datetime, timedelta
from typing import List, Dict, Any
from tkinter import messagebox
from dataclasses import dataclass

_INVISIBLE = re.compile(r'[\u200B-\u200D\uFEFF]')   # ZWSP, ZWNJ, ZWJ, BOM
def strip_invisibles(s: str) -> str:
    """Remove zero-width and BOM characters."""
    return _INVISIBLE.sub('', s)


class TestResults:
    def __init__ (self, patient_id: str, staircase: List[bool], left_eye: List[bool], responses : List[int], ISIs : List [int] ):
        self.patient_id = patient_id
        self.staircase = staircase
        self.left_eye = left_eye
        self.responses = responses
        self.ISIs = ISIs
        if ISIs:
            self.all_ISIs = np.unique(np.concatenate(ISIs))
        else:
            self.all_ISIs = np.array([], dtype=float)
        self.no_trials = len (ISIs)

def extract_from_txt (filename, patient_ID, flags : list[bool]):#flags is LEye,REye,staircase,random
    resp_list = []##response and ISI lists, each contain an array per trial
    isi_list = []
    staircase_list = []
    left_eye_list = []
    expiry = 7 #7 days to do it

    with open (filename, "r",encoding="utf-8") as file:
        text=strip_invisibles(file.read()).replace("â€‹", "")#For some reason the .txt file
        #was adding this as a suffix to the patient IDs
        print (text)

    blocks = re.split(r'-{83,}\s*\n', text)#boundaries per block

    id_pattern = re.compile(
    rf'Participant ID:\s*{re.escape(patient_ID)}\b',
    re.IGNORECASE)
    start_pattern = re.compile(r'VisualTDT Test Started Timestamp:\s*(.+)', re.IGNORECASE)
    type_pattern  = re.compile(r'Test Type:\s*(.+)', re.IGNORECASE)
    eye_pattern   = re.compile(r'Eye Tested:\s*(.+)', re.IGNORECASE)
    resp_pattern  = re.compile(r'Response:\s*([0-2](?:\s*,\s*[0-2])*)', re.IGNORECASE)
    tdt_pattern   = re.compile(r'TDTArray:\s*([0-9\.]+(?:\s*,\s*[0-9\.]+)*)', re.IGNORECASE)
    inclusion_pattern = re.compile(r'Include test in analysis:\s*(YES|NO)\b', re.IGNORECASE)

    # Compute cutoff for last 7 days
    now = datetime.now()
    cutoff = now - timedelta(days=expiry)

    # Prepare output structure

    patient_blocks = []
    for i, blk in enumerate(blocks):
        if id_pattern.search(blk):
            if 'ERROR IN PREVIOUS RESULTS' not in blocks[i+1]:
                patient_blocks.append(blk)
    #now have a list of blocks to look at
    for i, blk in enumerate(patient_blocks):
        m_start = start_pattern.search(blk)
        m_type  = type_pattern.search(blk)
        m_inclusion = inclusion_pattern.search(blk)

        if m_inclusion.group(1).strip() == "NO":
            continue
        if m_type.group(1).strip() == "Stepped Test":
            staircase = True
        else: staircase = False
        m_eye   = eye_pattern.search(blk)
        if m_eye.group(1).strip() == "Left Eye":
            left_eye = True
        else: left_eye = False
        m_resp  = resp_pattern.search(blk)
        m_tdt   = tdt_pattern.search(blk)


        if not all([m_start, m_type, m_eye, m_resp, m_tdt]):
            continue
        
        #Filter for flags, removing unwanted trials
        if (flags[0] == False and left_eye == True) or (flags[1] == False and left_eye == False):
            continue
        if (flags[2] == False and staircase == True) or (flags[3] == False and staircase == False):
            continue

        # Parse timestamp
        start_str = m_start.group(1).strip()
        """try:
            ts = datetime.strptime(start_str, '%A, %B %d, %Y %I:%M:%S %p')
        except ValueError:
            ts = datetime.strptime(start_str, '%d/%m/%Y %H:%M:%S')

        # Discard if older than 7 days
        if ts < cutoff:
            continue """     #removing (temporarily?)

        # Convert arrays
        resp_array = np.array([int(x) for x in m_resp.group(1).split(',')])

        isi_array  = np.array([float(x) for x in m_tdt.group(1).split(',')]) * 1000
        if staircase == False:
            isi_array.sort()
        
        resp_array[resp_array == 2] = 1   #2 signifies did not respond, which signifies the trial finished necause the TDT had been found. Can assume the rest of the ISIs are 1s as they are
        #greater than the TDT


        staircase_list.append(staircase)
        left_eye_list.append (left_eye)
        resp_list.append(resp_array)
        isi_list.append(isi_array)

    test_results = TestResults (patient_ID, staircase_list, left_eye_list, resp_list, isi_list)

    if test_results.no_trials == 0:
        messagebox.showerror("Error", f"No data found for Participant ID: {patient_ID} within the last {expiry} days")
        raise ValueError(f"No data found for Participant ID: {patient_ID} within the last {expiry} days")
    

    return test_results
