#!/usr/bin/env python3

#######################################################################

# Basic ideas:
#
# 1 - Grab all miRNA info.
# 2 - Read through UTRs, getting those from one gene at a time.
# 3 - Identify miRNA sites for this gene.
# 4 - Group overlapping miRNA sites in different species into a group
#

import sys
import os
import re

USAGE = ""
FILE_FORMATS = ""
GROUP_NUM = 0
MIR_FAM_ID = ""
LAST_UTR_ID = ""
OUTPUT_THIS_GENE_THIS_MIR = []
MIR_ID_2_SEED = {}
MIR_ID_SPECIES = {}
MIR_TYPE_2_MATCH = {}
SPECIES_START_END = {}
SPECIES_START_END_2_MATCH = {}
SPECIES_TO_UTR = {}
SPECIES_START_END_REMOVED = {}
SPECIES_START_END_2_MATCH_REMOVED = {}
GROUP_NUM_TO_SITE_TYPES = {}
GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST = {}
GROUP_TYPES_LIST_2_GROUP_TYPE = {}
SITE_TO_GROUP_NUM = {}
GROUP_NUM_TO_SPECIES = {}
GET_MATCH = {}
SITE_ID_2_SITE_TYPE = {}
SITE_ID_2_LENGTH = {}

GET_MATCH = {1: 1, 2: 1, 3: 1, 4: 0, 5: 0, 6: 0}

FIND_SITES_ALL_SPECIES = 1
REQUIRED_OVERLAP = 2
BEG_UTR_MASK_LENGTH = 0
VERBOSE = 1

def checkArguments():
    # Check for input and output file arguments
    # Print info if there are any problems

    USAGE = """
    Description: Search for predicted miRNA targets
                 using the modified TargetScanS algorithm.

    USAGE:
            ./{} miRNA_file UTR_file PredictedTargetsOutputFile

    Required input files:
            miRNA_file    => miRNA families by species
            UTR_file      => Aligned UTRs

    Output file:
            PredictedTargetsOutputFile    => Lists sites using alignment coordinates (MSA and UTR)
    For a description of input file formats, type
            ./{} -h

    Author: George Bell, Bioinformatics and Research Computing
    Translator: Mahdi Anvari, University of Tehran
    Version: 6.0
    Copyright (c) The Whitehead Institute of Biomedical Research
    """.format(os.path.basename(__file__), os.path.basename(__file__))

    FILE_FORMATS = """
    ** Required input files:

    1 - miRNA_file    => miRNA families by species

            contains three fields (tab-delimited):
                    a. miRNA family ID/name
                    b. seed region (7mer) for this miRNA
                    c. species ID in which this miRNA has been annotated
            ex:
            let-7/98        GAGGUAG 10090
            let-7/98        GAGGUAG 9606

            A miRNA family that is present in multiple species
            should be represented in multiple lines, one for each species.

    2 - UTR_file      => Aligned UTRs

            contains three fields (tab-delimited):
                    a. Gene/UTR ID or name
                    b. Species ID for this gene/UTR (must match ID in miRNA file)
                    c. Aligned UTR or gene (with gaps from alignment)
            ex:
            BMP8B   9606    GUCCACCCGCCCGGC
            BMP8B   9615    -GUG--CUGCCCACC

            A gene will typically be represented on multiple adjacent lines.
    """

    if len(sys.argv) > 1 and sys.argv[1] == "-h":
        print(USAGE)
        print(FILE_FORMATS)
        sys.exit(0)
    elif len(sys.argv) < 4:
        print(USAGE)
        sys.exit(0)
    elif not os.path.exists(sys.argv[1]):  # miRNA file not present
        print("\nI can't find the file {}\n".format(sys.argv[1]))
        print("which should contain the miRNA families by species.\n")
        sys.exit()
    elif not os.path.exists(sys.argv[2]):  # UTR file not present
        print("\nI can't find the file {}\n".format(sys.argv[2]))
        print("which should contain the Aligned UTRs.\n")
        sys.exit()

    miRNAfile = sys.argv[1]
    UTRfile = sys.argv[2]
    coordsFile = sys.argv[3]

    if os.path.exists(coordsFile):
        answer = input("Should I over-write {} [yes/no]? ".format(coordsFile))
        if answer.lower() not in ['y', 'yes']:
            sys.exit()

    return (miRNAfile, UTRfile, coordsFile)

# MIRNA_FILE, UTR_FILE, COORDS_FILE = checkArguments()

def makeSeedMatchRegex(seedMatch):
    # Turn a seed match region into a Python regular expression

    seedMatchLength = len(seedMatch)
    seedMatchPattern = ""

    for index, seedMatchNt in enumerate(seedMatch):
        if index < seedMatchLength - 1:
            seedMatchPattern += f"{seedMatchNt}-{{0,}}"
        else:
            seedMatchPattern += seedMatchNt

    return seedMatchPattern


def get_seeds(seedRegion, MIR_FAM_ID):
    global MIR_TYPE_2_MATCH

    # Get the 7mer-m8 seed match (exactly the reverse complement)
    seed2 = seedRegion
    rseed2 = seed2[::-1].translate(str.maketrans("AUCG", "UAGC"))

    # Get the 6mer seed match (7mer-m8 minus the first position)
    rseed6 = rseed2[1:]

    # Get the 6mer-1a seed match (6mer seed match minus the first position and add A at the end)
    rseed5 = rseed6 + "A"

    # Get the 7mer-1A seed match (6mer and add A at end)
    rseed1 = rseed6 + "A"

    # Get the 8mer-A1 seed match (7mer-m8 seed match and add A at end)
    rseed3 = rseed2 + 'A'

    # Get the 8mer-U1 seed match (7mer-m8 seed match and add U at end)
    rseed4 = rseed2 + 'U'

    # Make regex for searching by adding potential gaps between each pair of nts
    MIR_TYPE_2_MATCH[MIR_FAM_ID] = {
        1: makeSeedMatchRegex(rseed1),
        2: makeSeedMatchRegex(rseed2),
        3: makeSeedMatchRegex(rseed3),
        4: makeSeedMatchRegex(rseed4),
        5: makeSeedMatchRegex(rseed5),
        6: makeSeedMatchRegex(rseed6)
    }

def readMiRNAs():
    global MIR_ID_2_SEED, MIR_ID_SPECIES, MIRNA_FILE

    with open(MIRNA_FILE, 'r') as mir_family_data:
        for line in mir_family_data:
            # let-7/98	GAGGUAG	10090
            line = line.strip()
            line = line.replace('\r', '')  # For Windows and Mac

            # Public data format
            MIR_FAM_ID, mirSeedRegion, mirSpeciesID = line.split('\t')

            # Convert from RNA to DNA if needed
            mirSeedRegion = mirSeedRegion.replace('T', 'U', -1)

            MIR_ID_2_SEED[MIR_FAM_ID] = mirSeedRegion

            # Make sure we know which miRNA family is present in each species
            MIR_ID_SPECIES[f"{MIR_FAM_ID}::{mirSpeciesID}"] = 1

    # Get patterns for search for (one for each)
    for MIR_FAM_ID in sorted(MIR_ID_2_SEED.keys()):
        get_seeds(MIR_ID_2_SEED[MIR_FAM_ID], MIR_FAM_ID)

# readMiRNAs()

'''
with open(COORDS_FILE, "w") as coords:
    # Print output file header
    header = "a_Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\n"
    coords.write(header)
'''

def get_site_type_keys():

    global SITE_ID_2_SITE_TYPE, SITE_ID_2_LENGTH
    # Convert site type ID into name
    SITE_ID_2_SITE_TYPE = {
        1: "7mer-1a",
        2: "7mer-m8",
        3: "8mer-1a",
        4: "8mer-1u",
        5: "6mer-1a",
        6: "6mer"
    }

    SITE_ID_2_LENGTH = {
        1: 7,
        2: 7,
        3: 8,
        4: 8,
        5: 6,
        6: 6
    }

# get_site_type_keys()

def get_matches(MIR_FAM_ID, speciesID, matchType):
    global GROUP_NUM
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    alignment = SPECIES_TO_UTR[speciesID]
    match = MIR_TYPE_2_MATCH[MIR_FAM_ID][matchType]

    for match in re.finditer(match, alignment):
        start = match.start() + 1
        end = match.end()
        matched_sub_alignment = match.group()
        length = len(matched_sub_alignment)

        # Link site to site type
        SPECIES_START_END[f"{speciesID}::{start}::{end}"] = matchType

        # Link site to actual match (sequence alignment region, possibly including gaps)
        SPECIES_START_END_2_MATCH[f"{speciesID}::{start}::{end}"] = matched_sub_alignment

        # Backtrack far enough for matches with gaps
        match_start = match.start()
        match_end = match.end()
        backtrack_length = match_end - match_start
        match_start -= backtrack_length - 1

def drop_this_site(species, start, end):
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    key = f"{species}::{start}::{end}"

    if key in SPECIES_START_END:
        # Keep record of what we deleted
        SPECIES_START_END_REMOVED[key] = SPECIES_START_END[key]
        SPECIES_START_END_2_MATCH_REMOVED[key] = SPECIES_START_END_2_MATCH[key]

        del SPECIES_START_END[key]
        del SPECIES_START_END_2_MATCH[key]

        # print(f"We're deleting {species}::{start}::{end}")
        return 1
    elif key in SPECIES_START_END_REMOVED:
        # We already deleted this site
        # print(f"We already deleted {species}::{start}::{end}")
        return 1
    else:
        # This site doesn't exist
        # Was it already deleted, or is there a gap that requires adjustment???
        return 0

def get_subset_coords(alignment):
    global USAGE
    global FILE_FORMATS
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH

    match_pos = None
    start_plus_one_offset = None
    start_plus_two_offset = None
    end_minus_one_offset = None

    match = re.search(r'^[^-]-*([^-])', alignment)
    if match:
        match_pos = match.start(1) + 1
        start_plus_one_offset = match_pos - 1

    match = re.search(r'^[^-]-*[^-]-*([^-])', alignment)
    if match:
        match_pos = match.start(1) + 1
        start_plus_two_offset = match_pos - 1

    match = re.search(r'([^-])-*[^-]$', alignment)
    if match:
        end_minus_one_offset = len(match.group(1)) - 1

    return start_plus_one_offset, start_plus_two_offset, end_minus_one_offset

def find_remove_match_subsets():
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    # Remove shorter matches that are a subset of longer matches
    # ex: 8mer-1a includes 7mer-m8 and 7mer-1a matches

    species_start_end_copy = dict(SPECIES_START_END)
    #for site_this_utr_this_species in SPECIES_START_END:
    # Iterate over the copied dictionary
    for site_this_utr_this_species in species_start_end_copy:
        # Your code for processing each item goes here

        # Extracting information
        species, start, end = site_this_utr_this_species.split("::")
        start = int(start)
        end = int(end)
        gapped_match = SPECIES_START_END_2_MATCH.get(site_this_utr_this_species)

        # Calculating offsets based on gapped match
        if gapped_match and '-' not in gapped_match:
            start_plus_one = start + 1
            start_plus_two = start + 2
            end_minus_one = end - 1
        elif gapped_match:
            start_plus_one_offset, start_plus_two_offset, end_minus_one_offset = get_subset_coords(gapped_match)
            start_plus_one = start + start_plus_one_offset
            start_plus_two = start + start_plus_two_offset
            end_minus_one = end - end_minus_one_offset
            # print("Match with gaps:", gapped_match, "(", site_this_utr_this_species, ")", "[", start_plus_one, ",", start_plus_two, ",", end_minus_one, "]")

        drop_site = 0

        if SPECIES_START_END.get(site_this_utr_this_species) == 1:  # 7mer-1a
            # Drop 6mer with same start position and
            #      6mer-1a with same end position
            if 6 in GET_MATCH:
                drop_site = drop_this_site(species, start, end_minus_one)  # 6mer
            if 5 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_one, end)  # 6mer-1a
        elif SPECIES_START_END.get(site_this_utr_this_species) == 2:  # 7mer-m8
            # Drop 6mer with same end position
            if 6 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_one, end)  # 6mer
        elif SPECIES_START_END.get(site_this_utr_this_species) == 3:  # 8mer-1a
            # Drop 7mer-m8 with same starting position and
            #      7mer-1a with same ending position and
            #      6mer-1a with same ending position and
            #      6mer starting one position later
            if 2 in GET_MATCH:
                drop_site = drop_this_site(species, start, end_minus_one)  # 7mer-m8
            if 1 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_one, end)  # 7mer-1a
            if 5 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_two, end)  # 6mer-1a
            if 6 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_one, end_minus_one)  # 6mer
        elif SPECIES_START_END.get(site_this_utr_this_species) == 4:  # 8mer-1u
            # Drop 7mer-m8 with same starting position and
            #      6mer starting one position later
            if 2 in GET_MATCH:
                drop_site = drop_this_site(species, start, end_minus_one)  # 7mer-m8
            if 6 in GET_MATCH:
                drop_site = drop_this_site(species, start_plus_one, end_minus_one)  # 6mer

def get_utr_coords(align, end, utr_type):

    utr_beg = align[:int(end)].replace('-', '')
    start = len(utr_beg)
    end = start + SITE_ID_2_LENGTH[utr_type] - 1
    return start, end

def get_utr_coords2(align, end, type):
    utr_beg = align[:end].replace('-', '')

    start = len(utr_beg)
    end = start + SITE_ID_2_LENGTH[type] - 1

    return start, end

def make_list_non_redundant(list_str, sep):
    # Convert a list string into a non-redundant array
    values = str(list_str).split(sep)
    # Use a set to keep only unique values
    unique_list = list(set(values))
    # Sort the unique values
    unique_list.sort()
    return unique_list

def make_list_non_redundant2(list_str, sep):
    values = list_str.split(sep)
    unique_list = list(set(values))
    unique_list.sort(key=lambda x: int(x) if x.isdigit() else x)
    return unique_list

def group_this_pair(site1, site2):
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH

    if not SITE_TO_GROUP_NUM.get(site1) or not SITE_TO_GROUP_NUM.get(site2):
        if SITE_TO_GROUP_NUM.get(site1):  # Site 1 already part of a group
            # Set site2 to the same group as site1
            SITE_TO_GROUP_NUM[site2] = SITE_TO_GROUP_NUM[site1]
        elif SITE_TO_GROUP_NUM.get(site2):  # Site 2 already part of a group
            # Set site1 to the same group as site2
            SITE_TO_GROUP_NUM[site1] = SITE_TO_GROUP_NUM[site2]
        else:
            # Increment the group number
            GROUP_NUM += 1
            SITE_TO_GROUP_NUM[site1] = GROUP_NUM
            SITE_TO_GROUP_NUM[site2] = GROUP_NUM

def summarize_print_groups_this_gene_this_mirna():
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE
    global COORDS_FILE

    group_num_to_site_types_list = {}
    group_to_info = {}

    for group_num in sorted(GROUP_NUM_TO_SITE_TYPES.keys()):
        species_this_group = make_list_non_redundant(GROUP_NUM_TO_SPECIES[group_num], ";")

        if GROUP_NUM_TO_SITE_TYPES[group_num]:
            site_types = str(GROUP_NUM_TO_SITE_TYPES[group_num]).split(";")
            unique_site_types = make_list_non_redundant(GROUP_NUM_TO_SITE_TYPES[group_num], ";")
            unique_site_types_names_list = "+".join([SITE_ID_2_SITE_TYPE[int(site_type)] for site_type in unique_site_types])

        else:
            unique_site_types_names_list = ""

        group_num_to_site_types_list[group_num] = unique_site_types_names_list
        group_to_info[group_num] = " ".join(species_this_group)

    # Sort array by "group number" field
    OUTPUT_THIS_GENE_THIS_MIR.sort(key=lambda x: int(x.split("\t")[7]))

    for data_one_site_this_gene_this_mir in OUTPUT_THIS_GENE_THIS_MIR:
        f = data_one_site_this_gene_this_mir.split("\t")
        #print(f)
        group_num_this_site = int(f[7])
        site_type_this_site = int(f[8])
        f[8] = SITE_ID_2_SITE_TYPE[int(site_type_this_site)]
        group_type = group_num_to_site_types_list[group_num_this_site]
        data_one_site_this_gene_this_mir = "\t".join(f)

        if GROUP_TYPES_LIST_2_GROUP_TYPE.get(group_type):
            group_type = GROUP_TYPES_LIST_2_GROUP_TYPE[group_type]

        data_one_site_this_gene_this_mir += f"\t{group_type}\t{group_to_info[group_num_this_site]}"

        # Add the species list for this site type (but only if site type ne group type)
        if group_type != SITE_ID_2_SITE_TYPE[site_type_this_site]:
            nr_species_this_subset = make_list_non_redundant(GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[group_num_this_site].get(site_type_this_site, ""), " ")
            data_one_site_this_gene_this_mir += f"\t{' '.join(nr_species_this_subset)}"

        #print(data_one_site_this_gene_this_mir + "\n")
        with open(COORDS_FILE, "a") as coords:
            # Print output file header
            coords.write((data_one_site_this_gene_this_mir + "\n"))


def group_sites_this_gene_this_mirna():
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    for site1 in sorted(SPECIES_START_END.keys()):
        site1_species, site1_start, site1_end = site1.split("::")

        for site2 in sorted(SPECIES_START_END.keys()):
            site2_species, site2_start, site2_end = site2.split("::")

            if site1_species != site2_species:
                if site1_start == site2_start and site1_end == site2_end:
                    group_this_pair(site1, site2)
                elif site1_start == site2_start:
                    group_this_pair(site1, site2)
                elif site1_end == site2_end:
                    group_this_pair(site1, site2)
                elif site1_start > site2_start and site1_start <= site2_end:
                    num_overlap_nt = int(site2_end) - int(site1_start) + 1
                    if num_overlap_nt >= REQUIRED_OVERLAP:
                        group_this_pair(site1, site2)
                elif site1_end >= site2_start and site1_end < site2_end:
                    num_overlap_nt = int(site1_end) - int(site2_start) + 1
                    if num_overlap_nt >= REQUIRED_OVERLAP:
                        group_this_pair(site1, site2)
                elif (int(site1_start) > int(site2_start) and int(site1_end) < int(site2_end)) or \
                        (int(site2_start) > int(site1_start) and int(site2_end) < int(site1_end)):
                    group_this_pair(site1, site2)

    for this_site in sorted(SPECIES_START_END.keys()):
        annotated = ""
        site_all_info = this_site.split("::")
        species_this_site = site_all_info[0]

        if this_site not in SITE_TO_GROUP_NUM:
            GROUP_NUM += 1
            SITE_TO_GROUP_NUM[this_site] = GROUP_NUM

        if SITE_TO_GROUP_NUM[this_site] not in GROUP_NUM_TO_SITE_TYPES:
            GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]] = SPECIES_START_END[this_site]
        else:
            GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]] = str(GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]])
            GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]] += ";" + str(SPECIES_START_END[this_site])

        if SITE_TO_GROUP_NUM[this_site] not in GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST:
            GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]] = {SPECIES_START_END[this_site]: species_this_site}
        else:
            if SPECIES_START_END[this_site] not in GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][SPECIES_START_END[this_site]] = species_this_site
            else:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][SPECIES_START_END[this_site]] += " " + species_this_site

        if SPECIES_START_END[this_site] == "1":
            if GET_MATCH["6"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["6"] += " " + species_this_site
            if GET_MATCH["5"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["5"] += " " + species_this_site
        elif SPECIES_START_END[this_site] == "2":
            if GET_MATCH["6"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["6"] += " " + species_this_site
        elif SPECIES_START_END[this_site] == "3":
            if GET_MATCH["1"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["1"] += " " + species_this_site
            if GET_MATCH["2"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["2"] += " " + species_this_site
            if GET_MATCH["5"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["5"] += " " + species_this_site
            if GET_MATCH["6"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["6"] += " " + species_this_site
        elif SPECIES_START_END[this_site] == "4":
            if GET_MATCH["2"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["2"] += " " + species_this_site
            if GET_MATCH["6"]:
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]]["6"] += " " + species_this_site

        utr_start, utr_end = get_utr_coords(SPECIES_TO_UTR[site_all_info[0]], site_all_info[1], SPECIES_START_END[this_site])

        #GROUP_NUM_TO_SPECIES[SITE_TO_GROUP_NUM[this_site]] += species_this_site + ";"
        group_num = SITE_TO_GROUP_NUM.get(this_site)
        if group_num is not None:
            if group_num not in GROUP_NUM_TO_SPECIES:
                GROUP_NUM_TO_SPECIES[group_num] = ""
            GROUP_NUM_TO_SPECIES[group_num] += species_this_site + ";"
        else:
            print(f"Warning: {this_site} not found in SITE_TO_GROUP_NUM")


        if f"{MIR_FAM_ID}::{species_this_site}" not in MIR_ID_SPECIES:
            annotated = " "
        else:
            annotated = "x"

        OUTPUT_THIS_GENE_THIS_MIR.append(f"{LAST_UTR_ID}\t{MIR_FAM_ID}\t{species_this_site}\t{site_all_info[1]}\t{site_all_info[2]}\t{utr_start}\t{utr_end}\t{SITE_TO_GROUP_NUM[this_site]}\t{SPECIES_START_END[this_site]}\t{annotated}")

def group_sites_this_gene_this_mirna2():
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    # Can we drop this variable???  It doesn't seem to be used.
    pair_to_distance = {}

    ############### Check for position overlap between sites

    # Do an all vs. all comparison to identify overlaps
    for site1 in sorted(SPECIES_START_END.keys()):
        site1_species, site1_start, site1_end = site1.split("::")

        for site2 in sorted(SPECIES_START_END.keys()):
            site2_species, site2_start, site2_end = site2.split("::")

            # Skip comparison of same-species sites
            if site1_species != site2_species:
                ############  Choose combinations to give overlap  ############

                # Same start and end
                if site1_start == site2_start and site1_end == site2_end:
                    group_this_pair(site1, site2)
                    pair_to_distance[f"{site1} {site2}"] = int(site1_end) - int(site1_start) + 1
                # Same start
                elif site1_start == site2_start:
                    group_this_pair(site1, site2)
                    if int(site1_end) > int(site2_end):
                        pair_to_distance[f"{site1} {site2}"] = int(site1_end) - int(site1_start) + 1
                    else:
                        pair_to_distance[f"{site1} {site2}"] = int(site2_end) - int(site2_start) + 1
                # Same end
                elif site1_end == site2_end:
                    group_this_pair(site1, site2)
                    if int(site1_start) < int(site2_start):
                        pair_to_distance[f"{site1} {site2}"] = int(site1_end) - int(site1_start) + 1
                    else:
                        pair_to_distance[f"{site1} {site2}"] = int(site2_end) - int(site2_start) + 1
                # Offset one direction
                #     xxxxxxx
                #    xxxxxxx
                elif int(site1_start) > int(site2_start) and int(site1_start) <= int(site2_end):
                    num_overlap_nt = int(site2_end) - int(site1_start) + 1
                    if num_overlap_nt >= REQUIRED_OVERLAP:
                        group_this_pair(site1, site2)
                        pair_to_distance[f"{site1} {site2}"] = num_overlap_nt
                # Offset other direction
                #    xxxxxxx
                #     xxxxxxx
                elif int(site1_end) >= int(site2_start) and int(site1_end) < int(site2_end):
                    num_overlap_nt = int(site1_end) - int(site2_start) + 1
                    if num_overlap_nt >= REQUIRED_OVERLAP:
                        group_this_pair(site1, site2)
                        pair_to_distance[f"{site1} {site2}"] = num_overlap_nt
                # One within the other (with gaps)
                #      xxxxxxx          xxxxxxxxx
                #     xxxxxxxxx          xxxxxxx
                elif (int(site1_start) > int(site2_start) and int(site1_end) < int(site2_end)) or \
                     (int(site2_start) > int(site1_start) and int(site2_end) < int(site1_end)):
                    group_this_pair(site1, site2)
                    pair_to_distance[f"{site1} {site2}"] = num_overlap_nt

    for this_site in sorted(SPECIES_START_END.keys()):
        annotated = ""

        site_all_info = this_site.split("::")
        species_this_site = site_all_info[0]

        # This site is a group of 1, so no group info yet
        if not SITE_TO_GROUP_NUM.get(this_site):
            # If this group hasn't yet been assigned a number, give it one.
            GROUP_NUM += 1
            SITE_TO_GROUP_NUM[this_site] = GROUP_NUM

        if not GROUP_NUM_TO_SITE_TYPES.get(SITE_TO_GROUP_NUM[this_site]):
            # Start a list of site types for this group
            GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]] = SPECIES_START_END[this_site]
        else:
            # Add to the list of site types for this group
            GROUP_NUM_TO_SITE_TYPES[SITE_TO_GROUP_NUM[this_site]] += ";" + SPECIES_START_END[this_site]

        # Make a list of species in which a site type is found (in this group)
        GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][SPECIES_START_END[this_site]] += species_this_site + " "

        ###  If a wide site is present, its subset sites are also present
        if SPECIES_START_END[this_site] == 1:  # 7mer-1a
            if GET_MATCH.get(6):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][6] += species_this_site + " "  # 6mer
            if GET_MATCH.get(5):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][5] += species_this_site + " "  # 6mer-1a
        elif SPECIES_START_END[this_site] == 2:  # 7mer-m8
            if GET_MATCH.get(6):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][6] += species_this_site + " "  # 6mer
        elif SPECIES_START_END[this_site] == 3:  # 8mer-1a
            if GET_MATCH.get(1):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][1] += species_this_site + " "  # 7mer-1a
            if GET_MATCH.get(2):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][2] += species_this_site + " "  # 7mer-m8
            if GET_MATCH.get(5):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][5] += species_this_site + " "  # 6mer-1a
            if GET_MATCH.get(6):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][6] += species_this_site + " "  # 6mer
        elif SPECIES_START_END[this_site] == 4:  # 8mer-1u
            if GET_MATCH.get(2):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][2] += species_this_site + " "  # 7mer-m8
            if GET_MATCH.get(6):
                GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST[SITE_TO_GROUP_NUM[this_site]][6] += species_this_site + " "  # 6mer

        # Given the MSA coords, get the corresponding UTR coords
        utr_start, utr_end = get_utr_coords(SPECIES_TO_UTR[site_all_info[0]], int(site_all_info[1]), SPECIES_START_END[this_site])

        # Link each group to the species within it
        GROUP_NUM_TO_SPECIES[SITE_TO_GROUP_NUM[this_site]] += species_this_site + ";"

        # Is this miRNA annotated in this species?
        if not MIR_ID_SPECIES.get(f"{MIR_FAM_ID}::{species_this_site}"):
            annotated = " "
        else:
            annotated = "x"

        OUTPUT_THIS_GENE_THIS_MIR.append(f"{LAST_UTR_ID}\t{MIR_FAM_ID}\t{species_this_site}\t{site_all_info[1]}\t{site_all_info[2]}\t{utr_start}\t{utr_end}\t{SITE_TO_GROUP_NUM[this_site]}\t{SPECIES_START_END[this_site]}\t{annotated}")

def process_UTR_set():
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    # Look at each miRNA family
    for MIR_FAM_ID in sorted(MIR_ID_2_SEED.keys()):
        # Do one species' UTR at a time
        for speciesIDthisUTR in sorted(SPECIES_TO_UTR.keys()):
            # Is this miRNA in this species?
            # If so [or if we want to look anyway], look for sites

            if FIND_SITES_ALL_SPECIES or MIR_ID_SPECIES.get(f"{MIR_FAM_ID}::{speciesIDthisUTR}"):
                # print("Looking for", MIR_FAM_ID, "sites in species (", speciesIDthisUTR, ")....")

                for matchType in GET_MATCH.keys():
                    if GET_MATCH[matchType]:
                        get_matches(MIR_FAM_ID, speciesIDthisUTR, matchType)

                ###  Merge these types of sites when possible (new GB method)
                ###  Start by dropping 6mer sites that are included by 7mer sites
                ###  and then move to 7mer sites that are included in 8mer sites
                ###  site types: # 6mer 6mer-1a 7mer-1a 7mer-m8 8mer-1a 8mer-1u

                find_remove_match_subsets()

        # If there are any hit(s) for this miRNA in any species, group each orthologous set of sites
        num_sites_this_mirna_this_species = len(SPECIES_START_END)
        if num_sites_this_mirna_this_species:
            group_sites_this_gene_this_mirna()

            # Finish processing this gene/miR data
            summarize_print_groups_this_gene_this_mirna()

        # Empty out these UTR+MIR-specific variables after finishing this MIR_FAM_ID in this UTR
        OUTPUT_THIS_GENE_THIS_MIR.clear()
        SPECIES_START_END.clear()
        SPECIES_START_END_2_MATCH.clear()
        SPECIES_START_END_REMOVED.clear()
        SPECIES_START_END_2_MATCH_REMOVED.clear()

def read_utrs():
    global UTR_FILE
    global GROUP_NUM
    global MIR_FAM_ID
    global LAST_UTR_ID
    global OUTPUT_THIS_GENE_THIS_MIR
    global MIR_ID_2_SEED
    global MIR_ID_SPECIES
    global MIR_TYPE_2_MATCH
    global SPECIES_START_END
    global SPECIES_START_END_2_MATCH
    global SPECIES_TO_UTR
    global SPECIES_START_END_REMOVED
    global SPECIES_START_END_2_MATCH_REMOVED
    global GROUP_NUM_TO_SITE_TYPES
    global GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST
    global GROUP_TYPES_LIST_2_GROUP_TYPE
    global SITE_TO_GROUP_NUM
    global GROUP_NUM_TO_SPECIES
    global GET_MATCH
    global SITE_ID_2_SITE_TYPE
    global SITE_ID_2_LENGTH
    global FIND_SITES_ALL_SPECIES
    global REQUIRED_OVERLAP
    global BEG_UTR_MASK_LENGTH
    global VERBOSE

    with open(UTR_FILE, 'r') as utrs:
        for line in utrs:
            # NM_031304       hg18    -       GGCCCCAC
            line = line.strip()
            if line:  # Ignore empty lines
                # Public code format
                utr_id, this_species_id, this_utr = line.split('\t')

                if this_species_id:  # Skip consensus sequence (if present)
                    # Convert from RNA to DNA if needed
                    this_utr = this_utr.replace('T', 'U', -1)

                    # Mask beginning of UTR since miRNA can't target UTR right next to CDS
                    this_utr = 'N' * BEG_UTR_MASK_LENGTH + this_utr[BEG_UTR_MASK_LENGTH:]

                    if utr_id and utr_id != LAST_UTR_ID:
                        if VERBOSE:
                            # print("Processing", utr_id)
                            print("Processing", utr_id)

                        # Look for sites in this gene (UTR) and process the results
                        process_UTR_set()

                        # Empty out these UTR-specific variables after finishing this set of UTRs
                        SPECIES_TO_UTR.clear()
                        GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST.clear()
                        GROUP_NUM_TO_SITE_TYPES.clear()
                        GROUP_NUM_TO_SPECIES.clear()
                        SITE_TO_GROUP_NUM.clear()

                    # Add this UTR to the set
                    SPECIES_TO_UTR[this_species_id] = this_utr

                    LAST_UTR_ID = utr_id

    # Get the last one
    process_UTR_set()

# read_utrs()


MIRNA_FILE, UTR_FILE, COORDS_FILE = checkArguments()
readMiRNAs()


with open(COORDS_FILE, "w") as coords:
    # Print output file header
    header = "a_Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\n"
    coords.write(header)



get_site_type_keys()

read_utrs()
