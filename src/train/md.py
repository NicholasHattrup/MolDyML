import random
import numpy as np
import os
import mdtraj as md
from openbabel import openbabel as ob
from sklearn.cluster import KMeans


# A bunch of miscellanous methods for preparing md data for training 



def get_unique_frames(xyz_file, top_file=None, num_clusters=None):
    """
    Loads an XYZ file with multiple configurations into an mdtraj Trajectory object
    and returns the 30% most unique configurations based on clustering using the
    RMSD of the positions as the distance metric.

    Parameters
    ----------
    xyz_file : str
        The path to the XYZ file.

    Returns
    -------
    unique_frames : mdtraj.Trajectory
        The 30% most unique configurations as an mdtraj Trajectory object.

    """

    # Load the XYZ file
    if top_file is None:
        ob_conversion = ob.OBConversion()
        ob_mol = ob.OBMol()
        ob_conversion.ReadFile(ob_mol, xyz_file)

        # Convert the XYZ file to a PDB topology file
        top_file = f'{xyz_file[:-4]}.pdb'
        ob_conversion.WriteFile(ob_mol, top_file)





    # Load the XYZ file into an mdtraj Trajectory object
    #topology = md.load('file.xyz', format='xyz')
    traj = md.load(xyz_file, top=top_file)
    # Initialize the RMSD matrix
    num_frames,_,_ = traj.xyz.shape
    rmsd_matrix = np.zeros((num_frames, num_frames))

    # Compute the RMSD between each pair of frames
    print('Calculating RMSD')
    for i in range(num_frames):
        rmsd = md.rmsd(traj[i:], traj[i])
        rmsd_matrix[i, i:] = rmsd
        rmsd_matrix[:,i] = rmsd_matrix[i,:] # # by symmetry of the rmsd 
    
    # Perform k-means clustering
    if num_clusters is None:
        num_clusters = int(np.floor(num_frames*.1)) # Choose large cluster number


    # Perform k-means clustering
    kmeans = KMeans(n_clusters=num_clusters, n_init='auto')
    print('Starting K_means Clustering')
    kmeans.fit(rmsd_matrix)
    # Find the most representative frame in each cluster
    most_representative_frames = []
    most_represenative_frames_indexes = []
    for cluster_id in range(num_clusters):
        # Get the indices of the frames in the current cluster
        cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
        
        # Compute the average RMSD between each pair of frames in the cluster
        rmsd_values = rmsd_matrix[cluster_indices][:, cluster_indices]
        avg_rmsd_values = np.mean(rmsd_values, axis=0)
        
        # Find the frame in the cluster with the lowest average RMSD to all other frames
        most_representative_frame_index = cluster_indices[np.argmin(avg_rmsd_values)]
        most_representative_frame = traj[most_representative_frame_index]
        most_representative_frames.append(most_representative_frame)
        most_represenative_frames_indexes.append(most_representative_frame_index)

    # Print the most representative frames
    return most_representative_frames, most_represenative_frames_indexes









def combine_pos_force_xyz(pos_file, force_file, output_file=None):
    """
    Combines a position and force XYZ file into an extended XYZ file.

    Parameters
    ----------
    pos_file : str
        The name of the position XYZ file.
    force_file : str
        The name of the force XYZ file.
    output_file : str, optional
        The name of the output extended XYZ file. If not provided, defaults to 'combined.xyz'.

    Returns
    -------
    None

    """
    # Load the position and force data using NumPy
    pos_data = []
    with open(pos_file, 'r') as f:
        lines = f.readlines()

    # Get the number of atoms and configurations from the position data
    num_atoms=int(lines[0])
    num_configs=len(lines) // (num_atoms+2)
    atomic_symbols=[]
    descriptions=[]
    for i in range(2,num_atoms+2):
        atomic_symbols.append(lines[i].strip().split()[0])
    i=0
    while i < len(lines):
        if lines[i].strip().isdigit():
                i += 2
                continue
        else:
            # Append the atomic symbol and coordinates to the position data list
            pos_data.append([float(coord) for coord in lines[i].strip().split()[1:]])
            i += 1
    force_data = []
    with open(force_file, 'r') as f:
        lines = f.readlines()
    i=0
    while i < len(lines):
        if lines[i].strip().isdigit():
                i += 2
                continue
        else:
            # Append the atomic symbol and coordinates to the position data list
            force_data.append([float(force) for force in lines[i].strip().split()[1:]])
            i += 1
    

    pos_data=np.asarray(pos_data)
    force_data=np.asarray(force_data)


    # Check that the number of atoms and configurations is the same in both files
    if pos_data.shape != force_data.shape:
        print('Error: number of atoms or configurations does not match')
        exit()


    # Set the output file name
    if output_file is None:
        output_file = 'combined.xyz'

    # Combine the position and force data into an extended XYZ file
    with open(output_file, 'w') as f:
        # Loop over each configuration
        for config in range(num_configs):
            # Write the number of atoms and a comment line for this configuration
            f.write('{}\n'.format(num_atoms))
            f.write('Configuration {}\n'.format(config+1))

            # Loop over each atom and write its data to the file
            for atom in range(num_atoms):
                # Compute the starting index for this atom in the position and force data
                pos_start = atom+config*num_atoms
                force_start = atom+config*num_atoms

                # Write the atom symbol, position coordinates, and force coordinates
                f.write('{} {} {} {} {} {} {}\n'.format(atomic_symbols[atom],pos_data[pos_start, 0], pos_data[pos_start,1], 
                                                         pos_data[pos_start,2], force_data[force_start,0], force_data[force_start,1], 
                                                         force_data[force_start, 2]))




def concatenate_xyz_files(directories, substring, extension='.xyz', output_file=None):
    """
    Recursively searches a list of directories for XYZ files matching a substring,
    and concatenates them into one large XYZ file.

    Args:
        directories (str or list): A str or list of directory paths to search.
        substring (str): A substring to search for in the file names.
        extension (str): The file extension to search for (e.g. '.xyz').

    Returns:
        None.
    """

    if isinstance(directories, str):
        directories = [directories]

    # If output_file is not provided, create a default file name in the current directory
    if output_file is None:
        output_file = '{}{}'.format(substring, extension)

    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Loop over each directory in the list
        for directory in directories:
            # Recursively search for files matching the substring and extension in the directory
            for subdir, _, files in os.walk(directory):
                for file in files:
                    if file.endswith(extension) and substring in file:
                        # Read the XYZ file and append its contents to the output file
                        with open(os.path.join(subdir, file), 'r') as xyz:
                            f.write(xyz.read())






