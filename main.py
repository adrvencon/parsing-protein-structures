from pdb_parser import download_pdb, parse_pdb, report

if __name__ == "__main__":
    for pdb_id in ["1RBP", "1QRE", "1DGF"]:
        data = parse_pdb(download_pdb(pdb_id), pdb_id)
        report(data)