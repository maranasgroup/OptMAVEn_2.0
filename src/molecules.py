""" This modules manages molecules. """

from collections import defaultdict, OrderedDict
import itertools
import os
import tempfile
import warnings

from Bio.PDB import NeighborSearch, MMCIFParser, PDBParser, Residue, Selection, Structure
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import numpy as np

import charmm
from console import disp
import klaus
import maps
import standards
import user_input


_MAX_RANDOM_RESIDUE = 2**30 - 1


class Molecule(object):
    """ This class represents a molecule and provides functions for file input/output and analyzing/manipulating the structure.
    Note that a Molecule never stores its structure as an attribute: rather only a filename to a PDB file of the structure is stored.
    This is because Molecules are often stored as attributes of Experiments, which are in turn pickled.
    If the Molecule stored its structure object, then the (usually large) structure object would be pickled along with the Experiment, consuming excessive disk space.
    By storing only the file path, Molecules not only reduce drive usage but also store information only in portable formats (PDB files, instead of pickle files).
    """
    def __init__(self, name, file_, experiment, exclude_hetero=standards.HetAtmInclude):
        # Set the name of the molecule.
        if name != "".join(name.split()) or os.path.sep in name:
            raise ValueError("Molecule names may not contain whitespace or {}.".format(os.path.sep))
        self.name = name
        # Set the file path of the molecule.
        # All molecules must have a valid file path, but the path does not need to exist when the molecule is created.
        if not isinstance(file_, str):
            raise TypeError("Expected str, got {}.".format(type(file_)))
        self.file = file_
        # Set other attributes.
        self.exclude_hetero = exclude_hetero
        self.excluded_residue_ids = None
        self.experiment = experiment
        self.dims = standards.AtomCoordOrder  # atom coordinates are in order of x, y, z
        self.ndim = len(self.dims)  # molecules are 3-dimensional
    
    def get_structure(self, file_name=None):
        """ Load a structure from a file. """
        # The file doesn't need to be the filename associated with the molecule.
        # Molecules can be used to load structures from any PDB or MMCIF file.
        if file_name is None:
            file_name = self.file
        structure = get_parser(file_name).get_structure(self.name, file_name)
        # OptMAVEn currently can only handle structures with one model.
        models = Selection.unfold_entities(structure, "M")
        if len(models) != 1:
            raise NotImplementedError("Optmaven cannot currently handle structures with {} models.".format(len(models)))
        # Optionally, exclude heteroatoms from the structure. Only do this the first time this molecule loads a structure.
        if self.exclude_hetero in [standards.HetAtmAsk, standards.HetAtmExclude] and self.excluded_residue_ids is None:
            self.excluded_residue_ids = dict()
            # Loop through the chains.
            for chain in Selection.unfold_entities(models[0], "C"):
                # Identify any heteroatoms.
                hetero_residues = [residue for residue in chain if residue.get_id()[0] != " "]
                hetero_residue_ids = [residue.get_id() for residue in hetero_residues]
                if len(hetero_residue_ids) > 0:
                    # If there are heteroatoms:
                    if self.exclude_hetero == standards.HetAtmAsk:
                        # If the user can decide whether or not to exclude them, then ask.
                        disp("Excluding hetero-residues from chain {}.".format(chain.get_id()))
                        self.excluded_residue_ids[chain.get_id()] = set(user_input.select_from_list("Please select hetero-residues to exclude from {}: ".format(self.name), hetero_residue_ids, names=map(residue_code, hetero_residues)))
                    elif self.exclude_hetero == standards.HetAtmExclude:
                        # Otherwise, exclude heteroatoms automatically.
                        self.excluded_residue_ids[chain.get_id()] = hetero_residue_ids
                    else:
                        raise ValueError("Bad value for heteroatom exclusion: {}".format(self.exclude_hetero))
        elif self.exclude_hetero not in standards.HetAtmOptions:
            raise ValueError("Bad value for heteroatom exclusion: {}".format(self.exclude_hetero))
        if self.excluded_residue_ids is not None:
            # The previous step just identified residues to exclude.
            # This step actually removes them from the structure.
            for chain in Selection.unfold_entities(models[0], "C"):
                exclude = self.excluded_residue_ids.get(chain.get_id())
                if exclude is not None:
                    for res_id in exclude:
                        try:
                            chain.detach_child(res_id)
                        except KeyError:
                            pass
        return structure

    def get_model(self):
        return Selection.unfold_entities(self.get_structure(), "M")[0]
 
    def get_file_base(self):
        return os.path.basename(self.file)

    def get_file_split(self):
        return os.path.splitext(self.get_file_base())

    def get_file_name(self):
        return self.get_file_split()[0]

    def get_chains(self):
        return Selection.unfold_entities(self.get_model(), "C")

    def file_exists(self):
        return os.path.isfile(self.file)

    def get_renumbered_structure(self, first_residue=1):
        """ Renumber the structure, starting numbering at a given residue. """
        structure = self.get_structure()
        renumber_structure(structure, first_residue)
        return structure

    def write_pdb(self, structure, file_name=None, selector=None, preserve_atom_numbering=False):
        """ Write a PDB file of a structure to a file. """
        # This function is usually used to convert a structure to a PDB file that does not yet exist.
        # For example, a Molecule could be created with a filename that does not exist.
        # A structure from another file could be passed to this function in order to create that non-existant file.
        # However, this function can also create a file with a different name than that of the molecule's file.
        if file_name is None:
            file_name = self.file
        write_pdb(structure, file_name, selector, preserve_atom_numbering=preserve_atom_numbering)
   
    def write_renumbered_pdb(self, file_name, first_residue=1, selector=None):
        """ Perform renumbering and then write a PDB. """
        self.write_pdb(self.get_renumbered_structure(first_residue), file_name, selector)
        
    def relax(self, relaxed_file=None, in_place=True, ignore_solvation_initially=True):
        """ Perform a structural relaxation. """
        self.relax_CHARMM(relaxed_file, in_place, ignore_solvation_initially)
    
    def relax_CHARMM(self, relaxed_file=None, in_place=True, ignore_solvation_initially=True):
        """ Perform a structural relaxation with CHARMM.
        relaxed_file: the name of the file in which to store the relaxed structure.
        in_place: should the original file of the molecule be overwritten with the relaxed structure?
        The default is to relax in place, that is, to load the structure from the molecule's file, relax it, and then overwrite the original structure. """
        # In order for the relaxation to not be in place, a name for the relaxed file must be given.
        if relaxed_file is None:
            if not in_place:
                raise ValueError("A relaxation must be in place if no alternate relaxed file is given.") 
            relaxed_file = self.file
        # Perform the relaxation.
        with charmm.Relaxation(self.experiment, [self], relaxed_file, ignore_solvation_initially=ignore_solvation_initially) as relax:
            pass
        # If the relaxation was in place and a relaxed file was given, change the filename of the molecule to that of the relaxed file.
        # This way, all subsequent functions of this molecule will use that relaxed structure by default.
        if in_place:
            self.file = relaxed_file

    def disassemble_into_chains(self, chain_groups=None, prefix=None, structure=None, directory=None):
        """ Disassemble a molecule into its chains. Each chain becomes a Molecule. Return a list of these Molecules. """
        # The structure can be passed, or it can come from the Molecule's file.
        if structure is None:
            structure = self.get_structure()
        # The directory into which to save the file for each chain: by default, the same directory as the molecule's file.
        if directory is None:
            directory = os.path.dirname(self.file)
        # A prefix for the files of each chain: by default, the name of the molecule file (minus the extension).
        if prefix is None:
            prefix, suffix = os.path.splitext(os.path.basename(self.file))
        # OptMAVEn can currently handle only structures with one model.
        models = structure.get_list()
        if len(models) > 1:
            raise NotImplementedError("Optmaven cannot currently handle Molecules with {} models.".format(len(models)))
        # List the chains.
        chains = (models[0]).get_list()
        # Optionally, instead of fully disassembling the molecule, chains can be grouped together.
        # chain_groups is a list of lists of chain ids; the chains in each list will be grouped together in the output.
        # For example, chain_groups = [["A", "B"], ["H", "L"]] would create two molecules comprising chains A and B and chains H and L.
        if chain_groups is None:
            # If no chain groups are given, put each chain in its own group.
            chain_groups = [[chain.get_id()] for chain in chains]
        chain_molecules = list()
        # Loop through the chain groups.
        for chain_group in chain_groups:
            # Make a selector that will select only chains in the group.
            selector = SelectChains(chain_group)
            # Name the molecule by joining the chain ids.
            chain_ids = "".join(chain_group)
            chain_name = "{}_chain{}".format(self.name, chain_ids)
            file_name = os.path.join(directory, "{}_chain{}.pdb".format(prefix, chain_ids))
            # Write the PDB of only the chains in the group; restart numbering each molecule at residue 1.
            self.write_pdb(structure, file_name, selector, preserve_atom_numbering=False)
            # Create a Molecule for this chain group.
            chain_molecule = Molecule(chain_name, file_name, self.experiment)
            chain_molecules.append(chain_molecule)
        return chain_molecules

    def disassemble_for_CHARMM(self, structure=None, directory=None):
        """ Perform a special kind of disassembly for CHARMM:
        All chains go to their own molecule.
        Numbering is preserved when the files are written, i.e. each chain retains its original numbers, which CHARMM requires.
        """
        # Get the structure.
        if structure is None:
            structure = self.get_structure()
        # Get the directory into which the chain files will go.
        if directory is None:
            directory = self.directory
        # First, renumber the structure to ensure that atom and residue numbers start at 1 and that there are no numbering gaps.
        # Otherwise, CHARMM will crash.
        renumber_structure(structure)
        # Rename all of the residues so that they correspond to the CHARMM naming conventions.
        his_to_hsd(structure)
        # Ensure there is only one model.
        models = structure.get_list()
        if len(models) > 1:
            raise NotImplementedError("Optmaven cannot currently handle Molecules with {} models.".format(len(models)))
        # Loop through the chains in the model.
        chains = (models[0]).get_list()
        files = list()
        for chain in chains:
            # Make a selector to select only atoms in the chain.
            _id = chain.get_id()
            selector = SelectChains([_id])
            base_name = charmm.make_pdb_name(_id, "in")
            file_name = os.path.join(directory, base_name)
            # Write a PDB containing only atoms in the chain, preserving numbering.
            self.write_pdb(structure, file_name, selector, preserve_atom_numbering=True)
            files.append(file_name)
        # Return the chain objects (not Molecules, as in disassemble_into_chains) and a list of files.
        return chains, files

    def assemble_from_CHARMM(self, chain_ids, files, residue_numbers_by_model):
        """ After running CHARMM, assemble the output files into the molecule.
        This process is different from a simple merge because the molecules lose their chain ids and original residue numbering within CHARMM, so these ids must be added back. 
        """
        assembled_structure = None
        assembled_model = None
        # Loop through the chains to re-assemble.
        for _id, _file in zip(chain_ids, files):
            # Construct a structure from the chain's file.
            structure = self.get_structure(_file)
            # Get the single model from the structure.
            models = structure.get_list()
            if len(models) > 1:
                raise ValueError("A file output from CHARMM should have 1 model, not {}.".format(len(models)))
            model = models[0]
            # Get the original residue numbers in the model.
            residue_numbers_by_chain = residue_numbers_by_model[model.get_id()]
            # Get the single chain in the model.
            chains = model.get_list()
            if len(chains) > 1:
                raise ValueError("A file output from CHARMM should have 1 chain, not {}.".format(len(chains)))
            chain = chains[0]
            chain.id = _id
            # Get the original residue numbers in the chain.
            residue_numbers = residue_numbers_by_chain[_id]
            # Renumber the chain using the original residue numbers.
            renumber_chain(chain, residue_numbers=residue_numbers)
            # Add the chain to the assembled structure.
            if assembled_structure is None:
                # If no chain has been added yet, initialize the structure and model to the current structure object and its model.
                assembled_structure = structure
                assembled_model = assembled_structure.get_list()[0]
            else:
                # Otherwise, remove the current chain from its model and add it to the model that is being assembled.
                chain.detach_parent()
                assembled_model.add(chain)
        # Write a file of the assembled structure.
        self.write_pdb(assembled_structure)
        return assembled_structure

    def get_coords(self):
        """ Return the coordinates of the atoms in the molecule. """
        return self.get_atom_array()[self.dims].view(dtype=np.float).reshape(-1, self.ndim)

    def get_center(self):
        """ Return the geometric center of the molecule. """
        return np.mean(self.get_coords(), axis=0)

    def position_antigen(self, zAngle, x, y, z, file_name=None, in_place=False):
        """ Move an antigen molecule so that its epitope geometric center is at the given x, y, z coordinate and the antigen faces in the zAngle direction.
        Also rotate the antigen so as to minimize the z coordinate of the epitope geometric center.
        file_name and in_place work as they do in self.relax
        """
        # An output file must be given if the repositioned molecule is not supposed to overwrite the existing file.
        if file_name is None:
            if not in_place:
                raise ValueError("An output file must be given if the repositioning is not in place.")
            file_name = self.file
        # Reposition the antigen.
        with klaus.PositionAntigen(self.experiment, self.file, file_name, zAngle, x, y, z) as p:
            pass
        # If the repositioning was done in place, set the file name to the name of the repositioned file.
        if in_place:
            self.file = file_name

    def get_atoms(self):
        """ Return a list of all of the atom objects in the molecule. """
        structure = self.get_structure()
        chains = (structure.get_list()[0]).get_list()
        atoms = [a for c in chains for r in c for a in r]
        return atoms
        
    def get_atom_array(self):
        """ Return a numpy recordarray of the atoms in the molecule.
        The records are 'C' (chain id), 'R' (residue number), 'A' (atom name), and 'x', 'y', and 'z' (atomic coordinates).
        """
        raw_atoms = [(a.get_parent().get_parent().get_id(), str(a.get_parent().get_id()[1]), a.get_id(), a.coord[0], a.coord[1], a.coord[2]) for a in self.get_atoms()]
        c_format = "S{}".format(max([len(a[0]) for a in raw_atoms]))
        r_format = "S{}".format(max([len(a[1]) for a in raw_atoms]))
        a_format = "S{}".format(max([len(a[2]) for a in raw_atoms]))
        atom_format = [("C", c_format), ("R", r_format), ("A", a_format), (self.dims[0], float), (self.dims[1], float), (self.dims[2], float)]
        atom_array = np.array(raw_atoms, dtype=atom_format)
        return atom_array

    def interchain_residue_contacts(self, chain_ids_1, chain_ids_2, radius):
        """ Generate a list of residue contacts between two chains. """
        all_chains = {chain.get_id(): chain for chain in self.get_chains()}
        selected_chains = {chain_id: chain for chain_id, chain in all_chains.items() if chain_id in chain_ids_1 + chain_ids_2}
        atoms = [atom for chain_id, chain in selected_chains.items() for atom in Selection.unfold_entities(chain, "A")]
        residue_contacts = NeighborSearch(atoms).search_all(radius, "R")
        classified_contacts = defaultdict(list)
        for contact in residue_contacts:
            chain_1, chain_2 = [residue.get_parent().get_id() for residue in contact]
            if chain_1 in chain_ids_1 and chain_2 in chain_ids_2:
                classified_contacts[(chain_1, chain_2)].append({chain_1: contact[0], chain_2: contact[1]})
            elif chain_2 in chain_ids_1 and chain_1 in chain_ids_2:
                classified_contacts[(chain_2, chain_1)].append({chain_1: contact[0], chain_2: contact[1]})
        return classified_contacts


class AntibodyAntigenComplex(Molecule):
    """ This Molecule represents a complex of an antibody and an antigen. """
    def __init__(self, name, file_, experiment):
        Molecule.__init__(self, name, file_, experiment)
        # Identify the light and heavy chain loci.
        chain_ids = [chain.get_id() for chain in self.get_chains()]
        self.light_chain = None
        self.heavy_chain = None
        self.antigen_chains = list()
        # Ensure that there is exactly one heavy chain (id must be H) and one light chain (id must be K or L).
        # All other chains are assumed to be part of the antigen. There must be at least one antigen chain.
        for _id in chain_ids:
            if _id in standards.MapsLightChains:
                if self.light_chain is not None:
                    raise ValueError("An Antibody may not have multiple light chains.")
                self.light_chain = _id
            elif _id in standards.MapsHeavyChains:
                if self.heavy_chain is not None:
                    raise ValueError("An Antibody may not have multiple heavy chains.")
                self.heavy_chain = _id
            else:
                self.antigen_chains.append(_id)
        if self.light_chain is None or self.heavy_chain is None or len(self.antigen_chains) == 0:
            raise ValueError("An AntibodyAntigenComplex must have heavy, light, and antigen chains, not {}".format(", ".join(chain_ids)))
        self.antibody_chains = [self.heavy_chain, self.light_chain]


class ProtoAntibody(object):
    """ A proto-antibody stores a set of six MAPs parts and an antigen position. The select parts step creates ProtoAntibodies. """
    def __init__(self, parts, position, energy, gap_penalty):
        # Initialize attributes.
        self.parts = dict()
        self.position = position  # antigen position
        self.energy = energy  # interaction energy
        self.gap_penalty = gap_penalty
        self.light_chain = None
        self.coords = None
        self.cluster = None
        for cdr, number in parts.items():
            # Ensure that lambda and kappa chains are not mixed.
            chain = maps.get_chain(cdr)
            if chain in maps.light_chains:
                if self.light_chain is None:
                    self.light_chain = chain
                elif self.light_chain != chain:
                    raise ValueError("An antibody may not contain both lambda and kappa chains.")
            # Ensure that the part exists.
            maps.join(cdr, number)
            # Store the part.
            self.parts[cdr] = str(number)
        # Ensure that all CDRs are present (with no extras).
        all_cdrs = sorted(maps.get_cdr_names(self.light_chain))
        self_cdrs = sorted(self.parts.keys())
        if self_cdrs != all_cdrs:
            raise ValueError("Expected CDRs {}, got {}".format(",".join(all_cdrs), ",".join(self_cdrs)))
        # Ensure there are no clashes between MAPs parts.
        [maps.clashing(part1, part2) for part1, part2 in itertools.combinations(self.list_parts(), 2)]

    def list_parts(self):
        """ List the six MAPs parts in the ProtoAntibody. """
        return [maps.join(cdr, number) for cdr, number in self.parts.items()]

    def translate_cdr(self, cdr):
        """ Translate a CDR to its canonical name: K -> L is the only name change. """
        chain, region = maps.split_cdr(cdr)
        if chain in maps.light_chains:
            return maps.translate_chain(cdr, self.light_chain)
        else:
            return cdr

    def get_namesake_parts(self):
        """ Translate the names of all parts to their canonical names (K -> L). """
        parts = OrderedDict()
        for cdr in standards.MapsNamesakeCdrs:
            self_cdr = self.translate_cdr(cdr)
            parts[cdr] = maps.join(self_cdr, self.parts[self_cdr])
        return parts

    def get_labeled_position(self):
        """ Label the antigen position, e.g. {'zAngle': 120.0, 'x': 5.0, 'y': -2.5, 'z': 15.0} """
        return OrderedDict([(label, coord) for label, coord in zip(standards.PositionOrder, self.position)])
    
    def set_cluster(self, number):
        """ Assign a k-means cluster. """
        self.cluster = int(number)

    def get_coords(self):
        """ Return the coordinate vector used in k-means clustering.
        This vector concatenates the coordinates of the six MAPs parts, the x, y, and z coordinates, and sin(zAngle) and cos(zAngle). """
        if self.coords is None:
            # Only assemble this vector once: the first time it is needed.
            coords = list()
            # Get the standard order of the coordinates.
            order = standards.CoordOrder
            # Get the value of each of these coordinates.
            for coord in standards.CoordOrder:
                if coord in standards.PositionCoordOrder:
                    # The coordinate is a positional coordinate (zAngle, x, y, z)
                    if coord in standards.PositionOrder:
                        # The coordinate is a Cartesian coordinate (x, y, z)
                        value = self.position[standards.PositionOrder.index(coord)]
                    elif standards.AngleLabel in coord:
                        # If the coordinate is a rotational coordinate, then compute the sine or cosine, as appropriate.
                        if standards.SinLabel in coord:
                            _coord = coord.replace(standards.SinLabel, "")
                            fxn = np.sin
                        elif standards.CosLabel in coord:
                            _coord = coord.replace(standards.CosLabel, "")
                            fxn = np.cos
                        else:
                            _coord = coord
                            fxn = lambda x: x
                        value = fxn(standards.degrees_to_radians(self.position[standards.PositionOrder.index(_coord)]))
                    else:
                        value = None
                    if value is not None:
                        coords.append(value)
                        continue
                try:
                    # The coordinate is a CDR.
                    # Try getting the CDR and the dimension number from the coordinate label.
                    cdr, dim = standards.to_maps_coord(coord)
                except TypeError:
                    pass
                else:
                    # Add the coordinate to the coordinate vector.
                    cdr = self.translate_cdr(cdr)
                    part = maps.join(cdr, self.parts[cdr])
                    coords.append(maps.get_coordinates(part, self.gap_penalty)[dim])
                    continue
                # If the coordinate cannot be identified, raise a ValueError.
                raise ValueError("Bad coordinate: {}".format(coord))
            # Save the coordinates so that they don't need to be recomputed.
            self.coords = np.array(coords)
        return self.coords

    def to_molecule(self, name, _file, experiment):
        """ Create a Molecule object from the ProtoAntibody. """
        # First create the antibody heavy and light chains.
        directory = tempfile.mkdtemp(dir=experiment.get_temp(), prefix="to_mol_")
        ab_name = "ab_temp"
        ab_file = os.path.join(directory, "ab_temp.pdb")
        parts = [maps.join(cdr, number) for cdr, number in self.parts.items()]
        with klaus.CreateAntibody(experiment, parts, ab_file) as p:
            pass
        antibody = Molecule(ab_name, ab_file, experiment)
        # Then position the antigen.
        ag_name = "ag_temp"
        ag_file = os.path.join(directory, "ag_temp.pdb")
        x = self.position[standards.PositionOrder.index(standards.xLabel)]
        y = self.position[standards.PositionOrder.index(standards.yLabel)]
        z = self.position[standards.PositionOrder.index(standards.zLabel)]
        zAngle = self.position[standards.PositionOrder.index(standards.zAngleLabel)]
        antigen = Molecule(ag_name, experiment.epitope_zmin_file, experiment)
        antigen.position_antigen(zAngle, x, y, z, file_name=ag_file, in_place=True)
        # Merge the antigen and antibody and save the structure as a PDB..
        _complex = merge([antibody, antigen], return_molecule=True, merged_name=name, merged_file=_file, write_pdb=True)
        # Remove the directory in which the complex was assembled.
        experiment.safe_rmtree(directory)
        # Return the merged antibody-antigen complex.
        return _complex

    def __iter__(self):
        return iter(self.get_coords())
    
    # All comparisons (=, <, >) are based on the interaction energy.

    def __eq__(self, other):
        if isinstance(other, ProtoAntibody):
            return self.energy == other.energy
        else:
            return self.energy == other

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        if isinstance(other, ProtoAntibody):
            return self.energy < other.energy
        else:
            return self.energy < other

    def __gt__(self, other):
        if isinstance(other, ProtoAntibody):
            return self.energy > other.energy
        else:
            return self.energy > other

    def __le__(self, other):
        return self == other or self < other

    def __ge__(self, other):
        return self == other or self > other


class SelectChains(Select):
    """ SelectChains subclasses the Select class and enables the selection of specific chains while writing a PDB. """
    def __init__(self, chain_ids):
        Select.__init__(self)
        if not isinstance(chain_ids, (list, tuple, set)):
            raise TypeError("Chain IDs must be a list, tuple, or set, not {}.".format(type(chain_ids)))
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return int(chain.get_id() in self.chain_ids)


def write_pdb(structure, file_name, selector=None, preserve_atom_numbering=False):
    """ Write a PDB file from a given structure with a given file_name. Optionally, write specific atoms with selector and preserve atom numbering. """
    writer = PDBIO()
    writer.set_structure(structure)
    if selector is None:
        writer.save(file_name, preserve_atom_numbering=preserve_atom_numbering)
    else:
        writer.save(file_name, selector, preserve_atom_numbering=preserve_atom_numbering)


def his_to_hsd(structure):
    """ Rename histidines as 'HSD' for CHARMM. """
    rename_residues(structure, {"HIS": "HSD"})


def hsd_to_his(structure):
    """ Rename histidines as 'HIS' for everything but CHARMM. """
    rename_residues(structure, {"HSD": "HIS"})


def rename_residues(structure, names):
    """ Rename a specific set of residues in a structure. """
    # names is a dict of {old_name: new_name}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname
                try:
                    residue.resname = names[resname]
                except KeyError:
                    pass


def renumber_structure(structure, first_residue=1):
    """ Renumber the residues in a structure. """
    resnum = first_residue
    for model in structure:
        resnum = renumber_model(model, resnum) + 1


def renumber_model(model, first_residue=1):
    """ Renumber the residues in a model. """
    resnum = first_residue
    for chain in model:
        resnum = renumber_chain(chain, resnum) + 1
    return resnum


def renumber_chain(chain, first_residue=1, residue_numbers=None):
    """ Renumber the residues in a chain.
    first_residue: the number at which to start numbering the chain
    residue_numbers: a list of numbers with which to renumber the chain
    """
    residues = chain.get_list()
    resids = [residue.get_id() for residue in residues]
    # Ensure that there are no residues with the same number.
    resids_set = set(resids)
    if len(resids) != len(resids_set):
        raise ValueError("Repeat residue IDs detected.")
    # If a residue_numbers list is given, its length must match the number of residues in the chain.
    if residue_numbers is not None and len(residue_numbers) != len(resids):
        raise ValueError("Residue numbers has length {} but residues has length {}.".format(len(residue_numbers), len(resids)))
    # Loop through the residues.
    for index, residue in enumerate(chain):
        if residue_numbers is not None:
            # If a residue number was given, use it.
            resnum_new = residue_numbers[index]
        else:
            # Otherwise, renumber all residues in numerical order, starting at first_residue.
            resnum_new = index + first_residue
        # Make a new resid for the residue.
        resid_new = get_renumbered_resid(residue, resnum_new)
        if resid_new == resids[index]:
            # Skip this residue if its new number matches the old.
            continue
        while resid_new in resids_set:
            # Biopython will raise an exception if a residue is renumbered such that its number matches the number of any existing residue.
            # If it matches, then temporarily renumber the matching residue to something random.
            matching_residue_index = resids.index(resid_new)
            if matching_residue_index <= index:
                raise ValueError("Residue numbers are out of order.")
            resid_temp = get_renumbered_resid(residue, np.random.randint(0, _MAX_RANDOM_RESIDUE))
            while resid_temp in resids_set:
                resid_temp = get_renumbered_resid(residue, np.random.randint(0, _MAX_RANDOM_RESIDUE))
            renumber_residue(residues[matching_residue_index], resid_temp[1])
            resids_set.remove(resid_new)
            resids_set.add(resid_temp)
            resids[matching_residue_index] = resid_temp
        # Renumber the residue.
        renumber_residue(residue, resnum_new)
        resids_set.remove(resids[index])
        resids_set.add(resid_new)
        resids[index] = resid_new
    return resnum_new


def get_renumbered_resid(residue, number):
    """ Make a new id for a residue by changing its number but preserving the other parts of the id. """
    return (residue.get_id()[0], number, residue.get_id()[2])


def renumber_residue(residue, number, increment=False):
    """ Change the number of a residue by changing its id. """
    if increment:
        residue.id = get_renumbered_resid(residue, residue.get_id()[1] + number)
    else:
        residue.id = get_renumbered_resid(residue, number)


def merge(molecule_list, return_molecule=False, merged_name=None, merged_file=None, write_pdb=False):
    """ Merge all of the molecules in a list into one Molecule. """
    # First, ensure that all of the items in molecules are actually Molecules.
    if not all([isinstance(molecule, Molecule) for molecule in molecule_list]):
        raise TypeError("Expected a list of Molecule objects.")
    # Merge all of the chains into one model in one structure.
    merged_structure = None
    exp = None
    # Loop through the molecules.
    for molecule in molecule_list:
        # Ensure that all molecules are part of the same Experiment.
        if exp is None:
            exp = molecule.experiment
        elif molecule.experiment != exp:
            raise ValueError("Cannot merge two Molecules from different Experiments.")
        # Get the structure and model of the Molecule.
        molecule_structure = molecule.get_structure()
        molecule_models = list(molecule_structure)
        if merged_structure is None:
            # If the merged structure hasn't been initialized, set it to the structure of the first Molecule.
            merged_structure = molecule_structure
            merged_model = molecule_models[0]
        else:
            # Otherwise, move each chain from its original structure into the merged structure.
            for chain in molecule_models[0]:
                chain_id = chain.get_id()
                # Ensure that there are no duplicate chain names.
                if chain_id in [merged_chain.get_id() for merged_chain in merged_model]:
                    raise ValueError("Attempted to merge two Molecules that each have a chain {}.".format(chain_id))
                chain.detach_parent()
                merged_model.add(chain)
    if return_molecule:
        # Return a Molecule object if desired.
        if not all([isinstance(arg, str) for arg in [merged_name, merged_file]]):
            raise ValueError("A name and file must be given to make a Molecule.")
        merged_molecule = Molecule(merged_name, merged_file, exp)
        if write_pdb:
            # Also write a PDB if desired.
            merged_molecule.write_pdb(merged_structure)
        return merged_molecule
    else:
        # Otherwise, return the merged structure itself.
        return merged_structure


def residue_code(residue):
    """ Turn a Residue object into a code of its number and heteroatom flag. """
    if isinstance(residue, Residue.Residue):
        return "{}{}".format(residue.get_id()[1], residue.get_id()[2]).strip()
    else:
        raise TypeError("Expected a Residue, got {}.".format(type(residue)))


def get_parser(file_):
    """ Get a parser appropriate to the file format. """
    try:
        # Try to get the file extension to determine the file format.
        file_base, ext = os.path.splitext(file_)
    except ValueError:
        raise ValueError("Cannot obtain extension of file {}".format(file_))
    else:
        try:
            # Use a parser appropriate for the file format.
            return {
                ".pdb": PDBParser(),
                ".cif": MMCIFParser()
            }[ext]
        except KeyError:
            raise ValueError("Unknown molecular file format: {}".format(ext))
