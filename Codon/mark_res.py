#2015 Xiaoju Zhang
from pymol import cmd

def mark_res():
    """ mark_res <pdb_file> <highlight_residue>
    select the chain and residue that ICP dipeptide locate and highlight it.
    """
    #list of pdb id, chain id, and residues informations for ICP. the lower case pdbs are from blast,
    #upercase PDBs are from SGD pdb homologous search.
    pdb_input_lists = [('1fa0', 'A', '177-178', '117-147'),
                       ('1u5t', 'A', '42-43'),
                       ('1yke', 'A', '103-104', '43-73'),
                       ('2b1e', 'A', '294-295', '234-264'),
                       ('2PFV', 'A', '294-295', '234-264'),
                       ('1MNM', 'C', '177-178', '117-147'),
                       ('2ckz', 'D', '11-12'),
                       ('2pk9', 'D', '150-151+245-246', '90-120+185-215'),
                       ('2pmi', 'D', '150-151+245-246', '90-120+185-215'),
                       ('2xfv', 'A', '7-8'),
                       ('3esl', 'A', '191-192', '131-161'),
                       ('3n7n', 'E', '20-21'),
                       ('3t5v', 'E', '288-289', '228-258'),
                       ('4bh6', 'A', '365-366', '305-335')]
    #iterate pdbs
    for pdb_info in pdb_input_lists:
        pdb_id = pdb_info[0]
        chain_id = pdb_info[1]
        res_id = pdb_info[2]
        up_region = pdb_info[3] if len(pdb_info) == 4 else '0'
        pdb_file = pdb_id + ".pdb"
        chain_selection = "chain " + chain_id
        res_selection = "res " + res_id
        up_region_selection = "res " + up_region

        #command chains to process the pdb files.
        cmd.load(pdb_file)
        cmd.remove("solvent")
        cmd.color("gray80")
        cmd.show("cartoon")
        cmd.hide("line")
        cmd.select("monomer", chain_selection)
        cmd.color("cyan", "monomer")
        cmd.select("icp_dipep", chain_selection + " and " + res_selection)
        cmd.select("up_region", chain_selection + " and " + up_region_selection)
        cmd.select("none")
        cmd.color("red", "icp_dipep")
        cmd.color("yellow", "up_region")
        cmd.zoom(chain_selection)
        # cmd.zoom("all")
        cmd.png(pdb_id + ".png")
        cmd.save(pdb_id + ".pse")

        cmd.delete("all")
cmd.extend('mark_res', mark_res)
