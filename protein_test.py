import atom
import residue
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt

###############################################################################################
#TODOs:
#   make option to exclude/include backbone atoms(50% completion)
#   csv data logging
#   matplotlib/excel/google_sheets graphing methods
#   other data analysis methods[*priority*: benchmark to new analysis]
###############################################################################################

def get_dist(coord1, coord2):
    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)

def size_dist(cliques):
    # For situation where you are testing multiple clique generation methods
    freq_arr = [0,0,0,0,0,0,0]
    for i in cliques:
        freq_arr[len(i)] += 1
    return freq_arr

class Protein:
    def __init__(self, name, file_path):
        self.name = name
        self.file_path = file_path
        self.residues = {} # Dict attribute containing Residue objects. Made at reading pdb file
        self.atom_resid_map = {} # Dict attribute for atomid=>resid map. Made at reading pdb file
        self.mcsc_map = {} # Dict attribute for atomid=>mc/sc identifier. Made at reading pdb file 
        self.atom_cliques = None # list of cliques using atomic coordinates
        self.centroid_cliques = None # list of cliques use residue centroid coordinates
        self.centroid_clique_frequency = None # Is this size frequency?
        self.atom_clique_frequency = None # Is this size frequency?
        with open(self.file_path) as pdb_file:
            for line in pdb_file:
                if line[0:4] == "ATOM": # Start reading in data after this line in pdb file
                    res_id = int(line[22:26].strip(" ")) # ID number
                    atom_id = int(line[6:11].strip(" ")) # ID number
                    res_name = line[17:20].strip(" ") # Amino-acid 3-letter code
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    symbol = line[76:78].strip(" ") # N, C, O
                    atom_name = line[12:16].strip(" ") # CA, C, CB, etc
                    coords = (coordx, coordy, coordz)
                    atm = atom.Atom(symbol, atom_name, atom_id, coords) # create an atom object instance
                    if res_id in self.residues: # If res_id a key in self.residues, then add new atom info
                        self.residues[res_id].add_atom(atm)
                        self.atom_resid_map[atom_id] = res_id
                        self.mcsc_map[atom_id] = atm.mcsc # atm obj has attr mcsc, a str, "mc"=main-chain,"sc"=side-chain 
                    else: # else, reached a new res_id, create a new residue object instance
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm]) # Residue object instance
                        self.atom_resid_map[atom_id] = res_id
                        self.mcsc_map[atom_id] = atm.mcsc
        for i in self.residues:
            self.residues[i].update_COM()
            
    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues
    
    def get_atomid_res_map(self):
        return self.atom_resid_map
    
    def get_mcsc_map(self):
        return self.mcsc_map

    def get_cliques(self, clique_type, exclude_backbone=False, distance_cutoff=6):
        if clique_type == "centroid":
            if self.centroid_cliques is None:
                self.generate_cliques("centroid", exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
            return self.centroid_cliques
        elif clique_type == "atom":
            if self.atom_cliques is None:
                self.generate_cliques("atom", exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
                return self.atom_cliques
            else:
                return self.atom_cliques
        
        else: raise Exception("invalid clique type")
    
    def get_clique_frequency(self, clique_type):
        if clique_type == "centroid":
            if self.centroid_clique_frequency is None:
                self.freq_analysis("centroid")
            return self.centroid_clique_frequency
        elif clique_type == "atom":
            if self.atom_clique_frequency is None:
                self.freq_analysis("atom")
            return self.atom_clique_frequency
        else: raise Exception("invalid clique type")
            
    def generate_centroid_cliques(self, exclude_backbone=False, distance_cutoff=6):
        centroids = []
        centroid_res = {}
        for i in self.residues:
            centroids.append(self.residues[i].get_centroid())
            centroid_res[self.residues[i].get_centroid()] = self.residues[i]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)

    def generate_atom_cliques(self, exclude_backbone=False, distance_cutoff=6):
        coords = []
        res_atomid_map = {}
        mcsc_map = {}
        for i in self.residues:
            for j in self.residues[i].get_atoms():
                coords.append(j.get_coords())

        coords_array = np.array(coords)
        del coords
        tri = scipy.spatial.qhull.Delaunay(coords_array)
        edges = []
        for n in tri.simplices:
            resid_set = []
            mcsc_set = []
            for i in n:
                resid_set.append(self.atom_resid_map[i])
                mcsc_set.append(self.mcsc_map[i])
            unique_res_set = set(resid_set)
            unique_atm_type_set = set(mcsc_set)
            print("{} : {}".format(unique_res_set, unique_atm_type_set))
            return
        
            if len(unique_res_set) <= 1:
                next
            if len(unique_atm_type_set) <= 1:
                next
                
            edge = sorted([n[0], n[1]])
            edges.append((coords_res[edge[0]], coords_res[edge[1]]))
            edge = sorted([n[0], n[2]])
            edges.append((coords_res[edge[0]], coords_res[edge[2]]))
            edge = sorted([n[0], n[3]])
            edges.append((coords_res[edge[0]], coords_res[edge[3]]))
            edge = sorted([n[1], n[2]])
            edges.append((coords_res[edge[1]], coords_res[edge[2]]))
            edge = sorted([n[1], n[3]])
            edges.append((coords_res[edge[1]], coords_res[edge[3]]))
            edge = sorted([n[2], n[3]])
            edges.append((coords_res[edge[2]], coords_res[edge[3]]))
            
        graph = nx.Graph(edges)
        self.atom_cliques = list(nx.find_cliques(graph))
        temp_arr = []
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = coords_resid[tuple(list(coords_array[self.atom_cliques[i][j]]))]
            self.atom_cliques[i] = list(set(self.atom_cliques[i]))
            self.atom_cliques[i].sort()
            if self.atom_cliques[i] not in temp_arr:
                temp_arr.append(self.atom_cliques[i])
        self.atom_cliques = np.array(temp_arr)
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = self.residues[self.atom_cliques[i][j]]

    def generate_cliques(self, clique_type, exclude_backbone=False, distance_cutoff=6):
        if clique_type == "centroid":
            self.generate_centroid_cliques(exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
        elif clique_type == "atom":
            self.generate_atom_cliques(exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
        else: raise Exception("Invalid clique type")

    def getMaxMinDistance(self, coords):
        max_dist = 0
        min_dist = 10000
        if len(coords) == 1: return 0, 0
        for i in range(len(coords)-1):
            for j in range(i+1, len(coords)):
                dist = get_dist(coords[i], coords[j])
                max_dist = max([max_dist, dist])
                min_dist = min([min_dist, dist])
        return max_dist, min_dist

    def distance_analysis(self, clique_type):
        if clique_type == "centroid": 
            for i in self.centroid_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                print(clique_min, clique_max)
        elif clique_type == "atom":
            for i in self.atom_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                print(clique_min, clique_max)
        else: raise Exception("Invalid clique type")

    def freq_analysis(self, clique_type):
        freq_arr = [0,0,0,0,0,0,0]
        if clique_type == "centroid":
            for i in self.centroid_cliques:
                freq_arr[len(i)] += 1
            self.centroid_clique_frequency = freq_arr
        elif clique_type == "atom":
            for i in self.atom_cliques:
                freq_arr[len(i)] += 1
            self.atom_clique_frequency = freq_arr
        else: raise Exception("Invalid clique type")
        print(freq_arr)
        temp ='''objects = (1, 2, 3, 4, 5, 6)
        y_pos = np.arange(len(objects))
        values = freq_arr[1:]
        plt.bar(y_pos, values)
        plt.xticks(y_pos, objects)
        plt.ylabel("count")
        plt.title("freq of {}-based clique sizes for {}".format(clique_type, self.name))
        plt.show()'''
        

    
#protein = Protein("4quv", "C:\\alpha\\4quv.pdb")
#print(protein.centroid_cliques[0][1].get_centroid())

