import numpy as np
import os
import subprocess
import shutil

# spatial 
# enter the resolutions
# resolutions must be entered like : ----- coarse to finer -----
resolutions= np.array([8,16,32])
# change dt:
dt = 1e-3
dt_max = 1e-3
t_final = 1.e0 
q_x = np.zeros((np.size(resolutions),2))
# first column is 2nd IBM and other is staircase IBM
def save_init(init_file_path):
    init_open = open(init_file_path,"r")
    init_data_save = init_open.read()
    init_open.close()
    return init_data_save
def clean_init(init_file_path,old_data):
    open(init_file_path, "w").close()
    write_old(init_file_path,old_data)
    return
def write_old(init_file_path,old_data):
    f = open(init_file_path,"w")
    f.write(old_data)
    f.close()
    return

def run_cases(init_file_path:str,resolutions:np.array,dt:float,dt_max:float,t_final:float,old_init_data):
    # make a directory to save the output files
    os.makedirs(str(os.getcwd()+"//richard_out"),exist_ok=True)
    for i in range(np.size(resolutions)):
        os.makedirs(str(os.getcwd()+f"//richard_out//res_{resolutions[i]}"),exist_ok=True)
        #open the file
        init_file = open(init_file_path,"r+")
        init_variables = init_file.read()
        #edit the variables for IBM_G and IBM
        init_variables = init_variables.replace("g%nx = 50;",(f"g%nx = {resolutions[i]};"))
        init_variables = init_variables.replace("g%ny = 50;",f"g%ny = {resolutions[i]};")
        init_variables = init_variables.replace("g%nz = 20;",f"g%nz = {resolutions[i]};")
        
        init_variables = init_variables.replace(f"g%dt = 1e-4",f"g%dt = {dt};")
        init_variables = init_variables.replace(f"g%dt_max = 1e-3",f"g%dtmax = {dt_max};")
        init_variables = init_variables.replace(f"g%t_final = g%dt*10000",f"g%t_final = {t_final};")
        # we add seek cause pyth adds values to the end now it goes up and replaces them
        # write the new init vars. inside
        init_file.seek(0)          
        init_file.write(init_variables)
        init_file.truncate()       
        init_file.close()
        # compile the files
        print("compile has begun...")
        subprocess.run(["bash","./compile.sh"])
        print("compile has completed...")
        # first run the IBM_G case
        print("running the 2nd Order IBM case...")
        subprocess.run(["./build_ibm2nd/main"])
        shutil.move("mf_data.txt",str(os.getcwd())+f"//richard_out//res_{resolutions[i]}//mean_flow_2nd.txt")
        # now run the IBM case
        print("running the Staircase IBM case...")
        subprocess.run(["./build_ibm/main"])
        shutil.move("mf_data.txt",str(os.getcwd())+f"//richard_out//res_{resolutions[i]}//mean_flow_stair.txt")
        #reset init.f90
        clean_init(init_file_path,old_init_data)
        

    return

def find_value(file_name:str,keyword:str):
    f = open(file_name,"r")
    for line in f:
        if keyword in line:
            value = line.split(":")[-1].strip()
            return float(value)
def org_calc_values(resolutions:np.array):
    for i in range(np.size(resolutions)):
        wdr = str(os.getcwd()+f"/richard_out/res_{resolutions[i]}/")
        wdr_2nd = wdr+"mean_flow_2nd.txt"
        wdr_st = wdr+"mean_flow_stair.txt"  
        q_x[i,0] = find_value(wdr_2nd,"mean flow")
        q_x[i,1] = find_value(wdr_st,"mean flow")
    return 
        
def extrapolate():
    # richardson extrapolation formula is from:
    # Joel H. Ferziger. Computational Methods for Fluid Dynamics
    # Ch. 3.9 Eq. 3.53
    p_2nd = np.log((q_x[1,0]-q_x[0,0])/(q_x[2,0]-q_x[1,0]))/np.log(2.0)
    p_stair = np.log((q_x[1,1]-q_x[0,1])/(q_x[2,1]-q_x[1,1]))/np.log(2.0)
    return p_2nd,p_stair 

init_data_save = save_init("init.f90")
run_cases("init.f90",resolutions,dt,dt_max,t_final,init_data_save)
org_calc_values(resolutions)
p_2nd,p_stair = extrapolate()
print(p_2nd,p_stair)
clean_init("init.f90",init_data_save)