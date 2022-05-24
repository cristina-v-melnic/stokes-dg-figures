import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

______________________________________________________________________________
# Path to the folder of the experiment
# CHANGE HERE for every new refinement simulation results.
______________________________________________________________________________

experiment_folder = '../Experiments/Refinement_Experiments_2022-01-18_12-17-59'
______________________________________________________________________________

# Parameters used for plotting
# No need to change unless the plotting style is suboptimal.
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
        r"\usepackage[utf8x]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{cmbright}",
    ]),
})

params = {'font.size': 21,
          'legend.handlelength': 2,
          'figure.figsize': [10.4, 8.8],
          'lines.markersize': 20.0,
          'lines.linewidth': 4.5
          # 'lines.dashed_pattern': [3.7, 1.6]
          }

plt.rcParams.update(params)

_______________________________________________________________________________
# The section with all important functions starts here.

def check_directory(dir):
    """ Check if the directory exists and, if not, then create it. """
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        pass


def get_errors_from_files(outfiles_dir, outfiles, keys):
    '''
    For the outfiles in a specific directory returns a multidimensional array containing
    errors of the types speficied in keys by the string that has to be found in the outfiles, 
    for example keys = ["L2(u)", "L2(div(u))", "L2(p)"].
     	
    :param outfiles_dir: path (str) to the folder of a specific velocity space
    :param outfiles: list (str) of all files with outputs of parmoon calculations
    for a certain velocity space and example
    :param keys: list (str) of keywords of error types to look for in the files
    
    :return: an array (float) of shape [len(keys), lenth_error_list]
    '''
    # All errors from the outfiles in the outfiles_dir(ectory) are extracted in one list.
    errors = []

    for file in outfiles:
        with open(outfiles_dir + '/' + file, 'r') as f:
            found = 0
            text = f.readlines()
            
            # For each line check if the specific key is found by looping over them. 
            for line in text:
                for key in keys:
                    if key in line and not(key+"_boundary" in line):
                        found = 1

                        if "nan" in line:
                            print("A nan value found in {}".format(file))

                        for t in line.split():
                            try:
                                errors.append(float(t))
                            except ValueError:
                                pass
            
            # The condition below will tell us that a value of a specific
            # key is missing in a file, which means that something has
            # gone wrong in the simulation, so one must go to the script to check.
            if (found != 1):
                print("'{}' could not be found in {}".format(key, file))

    # Finding the nr of values for each type of error.
    lenth_error_list = int(len(errors)/len(keys))
    if (len(outfiles) != lenth_error_list):
        print("Attention! Not all files in the output folders have readable content.")

    # Reshaping the array to differentiate between types of errors (u, p, div u).
    errors = np.transpose(np.array(errors).reshape(
        (lenth_error_list, len(keys))))
    return errors


def get_h_from_files(outfiles_dir, outfiles):
    '''
    For the outfiles of H1 conforming elements in a specific directory 
    returns an array containing h(min, max) values.
     	
    :param outfiles_dir: path (str) to the folder of a specific velocity space
    :param outfiles: list (str) of all files with outputs of parmoon calculations
    for a certain velocity space and example

    :return: an array (float) with values of the grid width h.
    '''
    if "TH" or "MINI" in outfiles_dir:
        h = []
        for file in outfiles:
            with open(outfiles_dir + '/' + file, 'r') as f:
                found = 0
                text = f.readlines()
                for line in text:
                    if "h(min, max)" in line:
                        found = 1
                        h.append(float(line.split()[-1]))
                if (found != 1):
                    print("'{}' could not be found in {}".format(
                        "h(min, max)", file))
        return np.array(h)
    else:
        print("Warning: h can only be read from outfiles of H1 conforming elements")
        pass


def generate_plots(examples, velocity_spaces, outfiles, keys, save_in_dir="Results", conv_order=[]):
    '''
    Generates plots of convergence histories for a certain examples in the 'save_in_dir' folder of the current path.
    
    :param examples: list (str) of examples used for simulations
    :param velocity_spaces: list (str) of velocity spaces used in simulations
    :param keys: list (str) of string values that indicate which data to extract from the outfiles
    :param save_in_dir: (str) the path to the folder where to save the results, create if doesn't exist
    :param conv_order: list (int 0,1,2..) specify to add h^i to the plot, if not leave the list empty (default)
    
    :return: 0
    '''

    # Create the necessary folders if they don't exist.
    check_directory(save_in_dir.split('/')[0])
    check_directory(save_in_dir)
    velocity_spaces = sorted(velocity_spaces)

    for example in examples:
        #print(example)
        h_outfiles = sorted(os.listdir(
            experiment_folder + "/" + example + "/" + "MINI" + "/" + "output_files"))
        #print(h_outfiles)
        h = get_h_from_files(experiment_folder + "/" + example +
                             "/" + "MINI" + "/" + "output_files", h_outfiles)
        #print("h={}".format(h))
        h = h

        for i in range(len(keys)):
            for vel_spc in velocity_spaces:
                outfiles_dir = experiment_folder + '/' + \
                    example + '/' + vel_spc + '/' + 'output_files'
                outfiles = sorted(os.listdir(
                    experiment_folder + '/' + example + '/' + vel_spc + '/' + 'output_files'))
                #print(outfiles)

                errors = get_errors_from_files(outfiles_dir, outfiles, keys)
                #print(errors)
                if "BDM" in vel_spc:
                       plt.loglog(h, errors[i], label=vel_spc, marker='D', color ='tab:blue', alpha=0.7)
                elif "SV4" in vel_spc:
                       plt.loglog(h, errors[i], label=vel_spc, marker='H', color ='tab:pink', alpha=0.7)       
                elif "SV" in vel_spc:
                       plt.loglog(h, errors[i], label=vel_spc, marker='H', color ='tab:green', alpha=0.7)
                elif "RT0" in vel_spc:
                       plt.loglog(h, errors[i], label=vel_spc, marker='v', color ='tab:pink', alpha=0.7)           
                elif "RT" in vel_spc:
                       plt.loglog(h, errors[i], label=vel_spc, marker='v', color ='tab:red', alpha=0.7)       
                else:
                       plt.loglog(h, errors[i], label=vel_spc, marker='o', color = 'tab:orange', alpha=0.7)
            # if len(conv_order) != 0:
             #   for order in conv_order:

            # plt.title(example + ",   " + keys[i])

            if keys[i] == "L2(u)":
                plt.ylabel("$||\mathbf{u}-\mathbf{u}_{h}||_{L^{2}(\Omega)}$")

                if len(conv_order) != 0:
                    order = conv_order[0]
                    plt.loglog(h, np.power(h, order), linestyle='dashed',
                               label="$O(h^{})$".format(order), color = 'tab:purple')
            elif keys[i] == "L2(div(u))":
                plt.ylabel(
                    r"$||\nabla \cdot \mathbf{u}_{h}||_{L^{2}(\Omega)}$")
            elif keys[i] == "L2(p)":
                plt.ylabel("$||p-p_{h}||_{L^{2}(\Omega)}$")
                if len(conv_order) != 0:
                    order = conv_order[1]
                    plt.loglog(h, np.power(h, order), linestyle='dashed',
                               label="$O(h^{})$".format(order), color = 'tab:purple')

            plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
            plt.xlabel("$h$")

            plt.savefig(save_in_dir + "/" + example + "_" + keys[i] + ".png")
            plt.clf()


def get_res_from_folders(outfiles_dir="../Experiments/Optimal_sigma_2021-11-15_19-25-33/SinCos.h", nr_refinement=1):
    '''
    Obtains the values of the DG parameter sigma and the residual of the iterative solver of the linear system of equations.
    These are taken from the outfiles for a certain example in the Results folder of the current path, specified as parameters
    of the function.
    
    :param outfiles_dir: list (str) of all the outfiles in outfiles_dir(ectory)
   
    :return: a list of arrays - one with residual values and one with corresponding sigma values (array[double], array[double])
    '''
    folders = sorted(os.listdir(outfiles_dir))
    res = []
    sigma = []
    for folder in folders:
        file_name = os.listdir(outfiles_dir + '/' + folder + '/output_files')
        with open(outfiles_dir + '/' + folder + '/output_files/' + file_name[0], 'r') as f:
            found_r = 0
            found_s = 0
            text = f.readlines()
            for line in text:
                if "face_sigma_DG:" in line:
                    found_s = 1
                    sigma.append(float(line.split()[1]))
                if "residual:" in line:
                    found_r = 1
                    res.append(float(line.split()[6]))

            if (found_r != 1):
                print("'{}' could not be found in {}".format("residual:", file))
            if (found_s != 1):
                print("'{}' could not be found in {}".format(
                    "face_sigma_DG:", file))

    return np.array(res), np.array(sigma)


def plot_sigma_res_refinement():
    res_1 = get_res_from_folders(
        "../Experiments/Optimal_sigma_2021-11-15_19-25-33/SinCos.h", nr_refinement=1)
    res_2 = get_res_from_folders(
        "../Experiments/Optimal_sigma_2021-11-15_19-53-12/SinCos.h", nr_refinement=2)
    res_3 = get_res_from_folders(
        "../Experiments/Optimal_sigma_2021-11-15_20-02-30/SinCos.h", nr_refinement=3)
    plt.scatter(np.log(res_1[1]), np.log(res_1[0]), label="n_ref = 1")
    plt.scatter(np.log(res_2[1]), np.log(res_2[0]), label="n_ref = 2")
    plt.scatter(np.log(res_3[1]), np.log(res_3[0]), label="n_ref = 3")

    plt.xlabel("log(sigma)")
    plt.ylabel("log(res)")
    plt.legend()
    plt.savefig("sincos_BDM3.png")
    plt.show()


def plot_sigma_res_refinement_BDM1():
    res_1 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref1_it_2021-11-19_13-10-00/SinCos.h", nr_refinement=1)
    res_2 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref1_it_2021-11-19_13-10-25/SinCos.h", nr_refinement=2)
    res_3 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref1_it_2021-11-19_13-10-41/SinCos.h", nr_refinement=3)
    print(res_1)
    print(res_3)
    plt.scatter(np.log(res_1[1]), np.log(res_1[0]), label="n_ref = 1")
    plt.scatter(np.log(res_2[1]), np.log(res_2[0]), label="n_ref = 2")
    plt.scatter(np.log(res_3[1]), np.log(res_3[0]), label="n_ref = 3")

    plt.xlabel("log(sigma)")
    plt.ylabel("log(res)")
    plt.legend()
    plt.savefig("sincos_BDM1_res.png")
    plt.show()


def plot_sigma_res_refinement_BDM1_it():
    res_1 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref1_it_2021-11-16_11-05-04/SinCos.h", nr_refinement=1)
    res_2 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref2_it_2021-11-16_11-04-43/SinCos.h", nr_refinement=2)
    res_3 = get_res_from_folders(
        "../Experiments/Optimal_sigma_BDM1_ref3_it_2021-11-16_10-57-09/SinCos.h", nr_refinement=3)

    plt.scatter(np.log(res_1[1]), np.log(res_1[0]), label="n_ref = 1")
    plt.scatter(np.log(res_2[1]), np.log(res_2[0]), label="n_ref = 2")
    plt.scatter(np.log(res_3[1]), np.log(res_3[0]), label="n_ref = 3")

    plt.xlabel("log(sigma)")
    plt.ylabel("log(res)")
    plt.legend()
    plt.savefig("sincos_BDM1_it10.png")
    plt.show()


def get_sigma_cond_nr(path="../../software/matlab_program/bin/files/search_1/matlab_results"):
    '''
    Attention! Edit the directory with results from MATLAB before calling.


    '''
    sigma = []
    c1 = []
    c2 = []
    c3 = []
    with open(path, 'r') as f:
        text = f.readlines()
        for line in text:
            if not("sigma" in line):
                terms = line.split()
                sigma.append(float(terms[0]))
                c1.append(float(terms[1]))
                c2.append(float(terms[2]))
                c3.append(float(terms[3]))
    print(sigma)
    print(c1)
    print(c2)
    print(c3)

    plt.scatter(np.log(sigma), np.log(c1), label="n_ref = 1")
    plt.scatter(np.log(sigma), np.log(c2), label="n_ref = 2")
    plt.scatter(np.log(sigma), np.log(c3), label="n_ref = 3")

    plt.xlabel("log(sigma)")
    plt.ylabel("log(cond_nr)")
    plt.legend()
    plt.savefig("sincos_BDM1_cond.png")
    plt.show()


def get_conv_history():
    '''
    Attention! Edit the global directory before calling.


    '''
    # A list of all example names in the experiment_folder.
    examples = sorted(os.listdir(experiment_folder))
    #print(examples)

    # A list of names of velocity spaces (Important! all examples have the same folders inside
    # if that weren't true, you would have to loop over examples to get the velocity spaces).
    velocity_spaces = sorted(os.listdir(experiment_folder + '/' + examples[0]))
    #print(velocity_spaces)
    velocity_spaces_dir = experiment_folder + '/' + examples[0]
    print(velocity_spaces_dir)
    for vel_spc in velocity_spaces:
        # Also all velocity space folder have the same nr of output files with the same name
        outfiles = sorted(os.listdir(experiment_folder + '/' +
                                     examples[0] + '/' + vel_spc + '/' + 'output_files'))
        #print(outfiles)

    # Keywords to look for error values.
        keys = ["L2(u)", "L2(div(u))", "L2(p)"]

    # Generate convergence histories for multiple velocity spaces on one graph 
    # and create a folder for storage for each group of spaces.
    
    # for set in range(0,3):
    #    velocity_spaces_per_graph = velocity_spaces[set:len(velocity_spaces):3]
    #    generate_plots(examples, velocity_spaces_per_graph, keys, save_in_dir = "Results/" + str(velocity_spaces_per_graph))

        examples = ["SinCos.h", "polynomial_solution.h"]
        #examples = ["SinCos.h"]
        velocity_spaces_per_graph = ["BDM1", "RT1", "MINI", "RT0"]
        generate_plots(examples, velocity_spaces_per_graph, outfiles,
                       keys, save_in_dir="Results/" + str(velocity_spaces_per_graph), conv_order=[2, 1])
        velocity_spaces_per_graph = ["BDM2", "RT2", "TH2", "SV2"]
        generate_plots(examples, velocity_spaces_per_graph, outfiles,
                       keys, save_in_dir="Results/" + str(velocity_spaces_per_graph), conv_order=[3, 2])
        velocity_spaces_per_graph = ["BDM3", "TH3", "RT3", "SV3", "SV4"]
        generate_plots(examples, velocity_spaces_per_graph, outfiles,
                       keys, save_in_dir="Results/" + str(velocity_spaces_per_graph), conv_order=[4, 3])


def get_table(outfiles_dir='../Experiments/Optimal_sigma_ref3_2022-01-26_09-31-03/polynomial_solution.h'):
    '''
    Attention! Edit the directory before calling.


    '''
    out_folders = sorted(os.listdir(outfiles_dir))
    keys = ["symmetry_DG:", "face_sigma_DG:", "L2(u)", "L2(div(u))", "H1-semi(u)", "L2(p)", "H1-semi(p)"]
    order = []
    sigmas = []
    names = []

    k = -1
    i = -1
    j = 0
    print(len(out_folders))
    for folder in out_folders:
        v_space = folder.split('_')[0]
        if v_space == "TH2":
            break
        if v_space not in order:
            order.append(v_space)
            k += 1
        file_name = os.listdir(outfiles_dir + '/' + folder + '/output_files')
        table = []
        names = []
        with open(outfiles_dir + '/' + folder + '/output_files/' + file_name[0], 'r') as f:
            print(file_name[0])
            found_key = 0
            j = 0
            text = f.readlines()
            for line in text:
                for key in keys:
                    if key in line:
                        found_key = 1
                        if key == "face_sigma_DG:":
                            sigma = line.split()[-1]
                            # if sigma not in str(sigmas):
                            sigmas.append(float(sigma))
                            table.append("{:.1e}".format(float(sigma)))
                            names.append(key)

                            # i += 1
                            # j = 0
                            # else:
                            #	pass

                        elif key == "residual:":
                            table.append("{:.2e}".format(
                                float(line.split()[6])))
                            names.append(key)
                            j += 1
                        # elif key == "L2(u) :":
                        #    table[k, j] = float(line.split()[2])
                        #    j += 1
                        # elif key == "L2(div(u)) :":
                        #    table[k, j] = float(line.split()[2])
                        #    j += 1
                        # elif key == "L2(p) :":
                        #    table[k, j] = float(line.split()[2])
                        #    j += 1
                        else:
                            # print(line)
                            # print(f)
                            table.append("{:.2e}".format(
                                float(line.split()[-1])))
                            names.append(key)
                            j += 1
            print(names)
            print(table)
            line = ''

            
            for t in table:
                if t == table[-1]:
                    line = line + ' ' + t 
                else:
                    line = line + ' ' + t 
                
            if  '-1' in str(table[0]):    
                with open('table_{}.txt'.format(table[0]), 'a') as the_file:
                    the_file.write(line + "\n")
            elif '0' in str(table[0]):    
                with open('table_{}.txt'.format(table[0]), 'a') as the_file:
                    the_file.write(line + "\n")
            else:
                with open('table_{}.txt'.format(table[0]), 'a') as the_file:
                    the_file.write(line + "\n")       
                    
            if (found_key != 1):
                print("At least a term could not be found in {}".format(f))
    print(order)
    #print(sigmas)

    line = ''
    for n in names:
        if n == names[-1]:
            line = line + ' ' + n 
        else:
            line = line + ' ' + n 

    with open('header_table_{}.txt'.format(table[0]), 'a') as the_file:
        the_file.writelines("\n" + line)


def get_robustness_exp(dir_name='../Experiments/Reynolds_nr_Scaling_2021-11-23_09-33-58/SinCos.h',  keys=["L2(div(u))"], save_in_dir="Pressure_Robustness", nr_ref=4):
    folders = sorted(os.listdir(dir_name))
    print(folders)

    errors = []
    h = []
    err1 = []

    vel_spaces = ['']
    Re = []
    count = 0
    
    for folder in folders:
        name = folder.split('_')
        if name[0] != vel_spaces[count]:
            vel_spaces.append(name[0])
            count += 1
            #print(Re)
            Re = []
        Re.append(float(name[-1].split('=')[-1]))

        outfiles_dir = dir_name + '/' + folder + '/output_files'
        outfiles = sorted(os.listdir(outfiles_dir))

        if "TH" in name[0]:
            print("!!! {}".format(name[0]))	
            h.append(get_h_from_files(dir_name + "/" +
                                      folder + "/" + "output_files", outfiles))
            err1.append(get_errors_from_files(outfiles_dir, outfiles, keys))
            print(err1)
        else:
            pass
            
	
        errors.append(get_errors_from_files(outfiles_dir, outfiles, keys))
        #print(folder)
    #print("h = {}".format(h))
    #print('\n')
    print("err1 = {}".format(err1))
    print(err1[0])
    print("shape = {}".format(np.shape(err1)))

    #plt.loglog(h[0], errors[i, :], label=vel_spc, marker='o')
    print("h = {}".format(h[0]))
    print(np.reshape(err1[0], np.size(h[0])))
    
    print(np.shape(errors))
    print(len(Re))
    print(vel_spaces)	
    print(len(vel_spaces))
    print(np.shape(errors))
    #emergency code
    #for i in range(1,len(vel_spaces)+1):
     #   for j in range(len(Re)):
     #       h_plot = h[0]
     #       ers_plot = np.reshape(errors[j + len(Re)*(i-1)], np.size(h[0]))
     #       #plt.loglog(h_plot, ers_plot, label="$Re = {}$".format(Re[j]), marker='o', alpha = 0.7)
            
    for i in range(2,4):
        for j in range(len(Re)):
            h_plot = h[0]
            ers_plot = np.reshape(err1[j + len(Re)*(i-2)], np.size(h[0]))
            plt.loglog(h_plot, ers_plot, label="$Re = {}$".format(Re[j]), marker='o', alpha = 0.7)

        plt.ylabel(r"$||\nabla \cdot \mathbf{u}_{h}||_{L^{2}(\Omega)}$")
        plt.xlabel("$h$")
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)

        #plt.savefig(save_in_dir + "/" + vel_spaces[i] + ".png")
        plt.savefig(save_in_dir + "/TH_{}".format(i) + ".png")
        plt.clf()

        print(outfiles)


if __name__ == '__main__':
    # print(get_res_from_folders("../Experiments/Optimal_sigma_2021-11-15_19-25-33/SinCos.h", nr_refinement = 1))
    # print(get_res_from_folders("../Experiments/Optimal_sigma_2021-11-15_19-53-12/SinCos.h", nr_refinement = 2))

    # commented Nov 22
    # plot_sigma_res_refinement_BDM1()
    # get_sigma_cond_nr()

    #get_conv_history()
    #get_table(outfiles_dir='../Experiments/Optimal_sigma_ref1_2021-12-08_06-26-47/polynomial_solution.h')
    get_table()

    #get_robustness_exp(dir_name='../Experiments/Reynolds_nr_Scaling_2021-11-29_19-28-37/SinCos.h', save_in_dir="Pressure_Robustness_sin", nr_ref = 7)
    #get_robustness_exp(dir_name='../Experiments/Reynolds_nr_Scaling_2021-11-29_19-28-37/polynomial_solution.h', save_in_dir="Pressure_Robustness_pol", nr_ref = 7)
