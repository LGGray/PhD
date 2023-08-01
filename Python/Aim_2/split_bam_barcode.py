import subprocess
import sys

# Get input file and barcode file names from command line arguments
IN = sys.argv[1]
BARIN = sys.argv[2]

# Read in barcodes from file
print("Reading in barcodes ...")
BARCODES = {}
with open(BARIN) as f:
    for line in f:
        barcode = line.strip()
        BARCODES[barcode] = 1
print("done!")
keys = list(BARCODES.keys())
print("Total number of barcodes:", len(keys))

# Read in BAM file and split reads by barcode
print("Reading BAM file ...")
ALL = {}
n = 0
b = 0
with subprocess.Popen(["samtools", "view", IN], stdout=subprocess.PIPE) as proc:
    for line in proc.stdout:
        n += 1
        line = line.decode().strip()
        if "CB:Z:" in line:
            barcodeID = line.split("CB:Z:")[1].split("\t")[0]
            if barcodeID in BARCODES:
                b += 1
                if barcodeID not in ALL:
                    ALL[barcodeID] = []
                ALL[barcodeID].append(line)
print("done!")

barcode = 'TTTGTCATCTGGGCCA-1'
ALL[barcode]

# Write out SAM files for each barcode
print("Printing SAM files ...")
for barcodeID, lines in ALL.items():
    out = barcodeID + ".sam"
    with open(out, "w") as f:
        f.write("\n".join(lines))
print("done!")

print("Total number of lines read in:", n)
print("Total number of reads passing:", b)