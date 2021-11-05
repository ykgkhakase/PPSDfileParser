import enum
import gzip
import os
import re
import sys

class MolBlcokElements:
	MolName = "MolName"
	Software = "Software"
	NumOfAtoms = "NumOfAtoms"
	MolMiscNumbers = "MolMiscNumbers" #TODO
	MiscNumbers = "MiscNumbers" # TODO
	NumOfBonds = "NumOfBonds"
	Atoms = "Atoms"
	Bonds = "Bonds"
	Properties = "Properties"
	Version = "Version"
	CoordX = "CoordX"
	CoordY = "CoordY"
	CoordZ = "CoordZ"
	Element = "Element"
	Index1 = "Index1"
	Index2 = "Index2"
	BondType = "BondType"
	OrderOfOccurrence = "OrderOfOccurrence"
	ExitCode = "ExitCode"

class MolBlockException(Exception):
	pass

_RegexpMEnd = re.compile(r"^M\s+END")
_RegexpTag = re.compile(r"^>\s+<(.+?)>")
# as alias
MBE = MolBlcokElements

def _emptyMolBlock() -> dict:
	return {
		MBE.OrderOfOccurrence:0,
		MBE.ExitCode:1,
		MBE.MolName: "",
		MBE.Software : "",
		MBE.NumOfAtoms: 0,
		MBE.NumOfBonds: 0,
		MBE.Atoms: [],
		MBE.Bonds: [],
		MBE.Properties : [],
	}

def _stripAnyBreakLines(text: str) -> str:
	return text.replace('\r', '').replace('\n', '')

def _openFile(path):
    if os.path.splitext(path)[1] == '.gz':
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')
	
def _recordBlock(fp, numOfSkippedMol=0):
	count = 0
	allline = fp.readlines()
	fp.close()
	fidx = 0
	if numOfSkippedMol > 0:
		for idx, line in enumerate(allline):
			if line.startswith('$$$$'):
				count += 1
				if count == numOfSkippedMol:
					fidx = idx + 1 # In the below for loop, read just after '$$$$'
					break
	record = []

	for line in allline[fidx:]:
		record.append(line)
		if line.startswith('$$$$'):
			count += 1
			yield count, record
			record.clear()

def _parseOneBlock(molBlock: list):
	rt = _emptyMolBlock()
	props = []
	flagMEND = False
	for idx, text in enumerate(molBlock):
		if idx == 0:
			rt[MBE.MolName] = _stripAnyBreakLines(text)
		elif idx == 1:
			rt[MBE.Software] = _stripAnyBreakLines(text)
		elif idx == 3:
			text = _stripAnyBreakLines(text)
			if text[-5:] != 'V2000':
				raise MolBlockException('Only V2000 format is supported.')
			rt[MBE.NumOfAtoms] = int(text[0:3])
			rt[MBE.NumOfBonds] = int(text[3:6])
			k = text[6:]
			rt[MBE.MolMiscNumbers] = k[0:-1]
		elif idx > 3 and idx <= 3 + rt[MBE.NumOfAtoms]:
			rt[MBE.Atoms].append(_parseAtomBlock(text))
		elif idx > 3 + rt[MBE.NumOfAtoms] and idx <= 3 + rt[MBE.NumOfAtoms] + rt[MBE.NumOfBonds]:
			rt[MBE.Bonds].append(_parseBondBlock(text))
		elif text.startswith('$$$$'):
			break
		elif idx > 3 + rt[MBE.NumOfAtoms] + rt[MBE.NumOfBonds]:
			# The other "M" properties are discarded.
			if flagMEND == False: 
				if _RegexpMEnd.match(text):
					flagMEND = True
				continue
			props.append(text)
		
		rt[MBE.Properties] = _parsePropertyLists(props)

	return rt

def _parseTag(text: str):
	a = _RegexpTag.match(text)
	return a.group(1)

def _parsePropertyLists(props: list) -> list:
	rt = []
	indexes = [idx for idx in range(len(props)) if props[idx].startswith('>')]
	for idx in range(0, len(indexes)):
		a1 = indexes[idx]
		if idx+1 == len(indexes):
			a2 = len(props)
		else:
			a2 = indexes[idx+1] - 1 
		rt.append({
			_parseTag(props[a1]) : [_stripAnyBreakLines(e) for e in props[a1+1:a2]]
		})
	return rt

def _parseAtomBlock(text: str) -> dict:
	k = _stripAnyBreakLines(text).split()

	return {
		MBE.CoordX : float(k[0]),
		MBE.CoordY : float(k[1]),
		MBE.CoordZ : float(k[2]),
		MBE.Element : k[3],
		MBE.MiscNumbers: [int(e) for e in k[4:]], 
	}

def _parseBondBlock(text: str) -> dict:
	text = _stripAnyBreakLines(text)
	return {
		MBE.Index1 : int(text[0:3]),
		MBE.Index2 : int(text[3:6]),
		MBE.BondType : text[6:9],
		MBE.MiscNumbers: [int(e) for e in text[9:].split()], 
	}
	pass

def SDFileParser(sdfile_path: str, maxNumOfMol = 0, numOfSkippedMol = 0):
	"""[Pure-Python implementation of Parser for SD file ]

	Args:
	    sdfile_path (str): [Filepaht of .sdf/.sdf.gz]
	    maxNumOfMol (int, optional): [maximum limit number of molecules]. If 0 is specifed, all molecules are read.
		numOfSkippedMol (int, optional): [skip reading N molecules from the head] If 0 is specifed,  any molecules are not skipped.

	Returns:
		list of molecule which consists of list/dict.
	"""
	rt = []
	fp = _openFile(sdfile_path)

	count = 0
	for idx, mb in _recordBlock(fp, numOfSkippedMol=numOfSkippedMol):
		try:
			r = _parseOneBlock(mb)
			r[MBE.ExitCode] = 0
		except ValueError:
			print(f"Value Error in #{idx} block of molecule.", file=sys.stderr)
			r[MBE.ExitCode] = 1
		count += 1
		r[MBE.OrderOfOccurrence] = idx
		rt.append(r)
		if maxNumOfMol != 0 and count == maxNumOfMol:
			break
	
	return rt

if __name__ == '__main__':

	def callSDFileParser(args):
		mols = SDFileParser(args.sdfile, 
			maxNumOfMol = args.maxNumOfMol,
			numOfSkippedMol= args.numOfSkippedMol)
		print(json.dumps(mols, indent=2), file=args.destFile)

	import json
	import argparse

	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()
	lp = subparsers.add_parser('tojson', help='show contents as JSON format')
	lp.add_argument('sdfile', action='store', help='a file path of ".sdf/.sdf.gz"')
	lp.add_argument('--maxNumOfMol', metavar='integer', action='store', type=int, default=0, help='Maximum number of molecules for output data')
	lp.add_argument('--numOfSkippedMol', metavar='integer', action='store', type=int, default=0, help='Number of molecules for skipping reading from the head record.')
	lp.add_argument('--destFile', action='store', default=sys.stdout, type=argparse.FileType('w', encoding='utf-8'))

	lp.set_defaults(func=callSDFileParser) 

	args = parser.parse_args()
	has_func = hasattr(args, 'func')
	if not has_func:
		parser.error('too few arguments')
		sys.exit(1)

	args.func(args)



