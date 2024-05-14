import pandas as pd
import numpy as np
import os
print('Enter the drive')
drive = input()
print('Enter the folder name:')
folder= input()
os.mkdir(''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\')
df = pd.read_excel('Parametric_X_Axial_Balanced_Shell.xlsx')
column=['filename','alpha','beta','gamma','tau','C_D','e_f_c','e_f_b','num_B_T','num_C_T']
df=df.reindex(column, axis='columns')
df.reset_index(inplace = True)
#Materials Input with Dimensions
import os
import math
E = 210000 #Modulus of elasticity
v = 0.3 #Poisson's ratio
for count in range(len(df)):
    inp_dir =count
    path = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'+str(inp_dir)+'/'
    directory = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'
    os.mkdir(path)
    alpha = df['alpha'][count]
    beta =df['beta'][count]
    gamma = df['gamma'][count]
    tau = df['tau'][count]
    C_D = df['C_D'][count]
    e_f_c = df['e_f_c'][count]
    e_f_b = df['e_f_b'][count]
    num_B_T = df['num_B_T'][count]
    num_C_T = df['num_C_T'][count]
    C_L = 0.5*alpha*C_D
    B_D = C_D*beta
    C_T = 0.5*C_D/gamma
    B_T = tau*C_T
    a_r_c = (1.64-(beta-0.2)*(1.64-0.75)*2)
    a_r_b = (0.8+(beta-0.2)*(1.2-0.8)*2)
    a_r_cc = 1/beta
    a_r_cs = (beta+(1.1-(beta-0.2)*0.3)*0.5*(1+beta))/(0.5*beta*3.14)
    if alpha<= 12:
        ratio_br = 1
        ratio_ch = 1
    else:
        ratio_br = alpha/8
        ratio_ch = alpha/8 
    ratio_br = 2
    ratio_ch = 2
    C_M_L =C_L/2
    C_R = C_D/2
    C_I_R = C_R - C_T
    B_R = B_D/2
    B_I_R = B_R - B_T
    B_M_R = B_R - 0.5*B_T
    C_M_R = C_R - 0.5*C_T
    B_L = round(C_M_L) 
    z = (1.1-(beta-0.2)*0.3)*0.5*(B_R+C_R)
    bcl = min(e_f_b*C_R,0.95*(B_L-C_R))
    z_1 = (e_f_c-z/C_R)*C_R
    theta_cut = 90-math.atan(C_R/(B_R+z-0.5*B_R))*180/math.pi
    ao = 0.2*math.sqrt(B_R*B_T)
    lrmin = max(4,0.4*C_T )
    if ao<=lrmin:
        a=lrmin
    else:
        a=ao
    lrmin = max(4,0.4*B_T )
    if ao<=lrmin:
        a_b=lrmin
    else:
        a_b=ao
    b_b = 0.65*math.sqrt(B_R*B_T)
    b_cc = 0.4*math.sqrt(math.sqrt(B_R*C_R*B_T*C_T))
    b_cs = math.pi*C_R/36
    e_s_c = C_T/num_C_T
    e_s_b = B_T/num_B_T
    num_cir_crown = round(0.5*math.pi*C_R/(e_s_c*a_r_cc))
    num_cir_saddle = round((B_R+z)/(e_s_c*a_r_cs))
    c= 1+round(a_r_c*(b_cc-a)/e_s_c)
    d= 1+round(a_r_b*(b_b-a_b)/e_s_b)
    num_weld_c = 1+round(a_r_c*a/(e_s_c))
    num_weld_b = 1+round(a_r_b*a_b/e_s_b)
    num_br_cut = round((2*d+num_weld_b)+a_r_b*(bcl+a_b-2*b_b)/e_s_b)
    num_ch_cut = round((2*c+num_weld_c)+a_r_c*((z+a-2*b_cs)/e_s_c))
    num_ch_cut_1 = round((a_r_c*((z_1)/e_s_c)))
    e_s_b_cut = (bcl+a-2*b_cs)/(num_br_cut-(num_weld_b+2*d)) 
    e_s_c_cut =  (z+a-2*b_cs)/(num_ch_cut-(num_weld_c+2*c))
    num_brace = round((B_L-C_R-bcl)/(ratio_br*e_s_b))
    num_chord = round((C_M_L-B_R-z-z_1)/(ratio_ch*e_s_c))
    e=round(2*B_I_R/e_s_c)
    num_ch_bt = round((1-(beta-0.2))*num_cir_crown)
    crown_r = round(B_R+(0.5*B_T),2)  
    weld_saddle = round((1-(math.asin(B_R/C_R)*180/math.pi)/45)*0.5*B_T ,1)
    weld_sad1 = max(0.05*B_T,weld_saddle)
    saddle_r = round(B_R + weld_sad1,2)
    dihedral_ang = 90+math.asin(B_R/C_R)*180/math.pi   
    jt_ang_cr= 46.5
    jt_ang_sd= min(54,37.5+(0.225*(dihedral_ang-50))) 
    jt_ang_sd= 37.5+(0.225*(dihedral_ang-50)) 
    with open(path +'model.txt','a+') as f:
      f.truncate(0)
      f.seek(0)
      f.write('/PREP7\n')
      f.write('BOPTN,NWARN,1\n') 
      f.write('ET,1,SHELL281 \n')
      f.write('sect,1,shell,, \n')
      f.write('secdata,'+str(C_T)+',1,0,3'+'\n')
      f.write('secoffset,MID \n')
      f.write('seccontrol,,,, , , , \n')  
      f.write('sect,2,shell,, \n')
      f.write('secdata,'+str(B_T)+',1,0,3'+'\n')
      f.write('secoffset,MID \n')
      f.write('seccontrol,,,, , , , \n') 
      f.write('MPTEMP,,,,,,,,\n')
      f.write('MPTEMP,1,0 \n')
      f.write('MPDATA,EX,1,,'+str(E)+'\n')
      f.write('MPDATA,PRXY,1,,'+str(v)+'\n')
      f.write('k,1,0,0,0'+'\n')
      f.write('k,2,0,0,'+str(C_M_L)+'\n')
      f.write('k,3,0,'+str(B_L)+',0'+'\n')
      f.write('l,1,2'+'\n')  
      f.write('l,1,3'+'\n') 
      f.write('circle,1,'+str(C_M_R)+',2'+'\n')  
      f.write('circle,1,'+str(B_M_R)+',3'+'\n')  
      f.write('ldele,4,9'+'\n')
      f.write('adrag,3,,,,,,,,,,,1'+'\n')
      f.write('adrag,10,,,,,,,,,,,2'+'\n')
      f.write('asba,2,1'+'\n')
      f.write('adele,4'+'\n')
      f.write('circle,1,'+str(C_M_R)+',2'+'\n')
      f.write('ldele,4,5'+'\n')  
      f.write('adrag,3,6,,,,,,,,,1,1'+'\n')
      f.write('lsla,u'+'\n') 
      f.write('ldele,all'+'\n') 
      f.write('allsel,all'+'\n') 
      f.write('aptn,all'+'\n')
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,-90 \n')
      f.write('wpoff,,,'+str(C_R+bcl)+'\n') 
      f.write('asbw,all'+'\n')
      f.write('wprota,,90 \n') 
      f.write('wpoff,,,'+str(B_R+z)+'\n')  
      f.write('asbw,all'+'\n') 
      f.write('allsel,all'+'\n')
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,-90 \n')  
      f.write('wpoff,,,'+str(C_R+a_b)+'\n')
      f.write('lsbw,17'+'\n')
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,14'+'\n')
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,16'+'\n')
      f.write('KWPAVE,17'+'\n')
      f.write('wpoff,,,'+str(a_b+0.5*C_T)+'\n')
      f.write('lsbw,18'+'\n') 
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,26'+'\n') 
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,18'+'\n') 
      f.write('WPCSYS,-1,0 \n')
      f.write('wpoff,,,'+str(B_R+a)+'\n')
      f.write('lsbw,24'+'\n')
      f.write('wpoff,,,'+str(b_cc-a)+'\n')
      f.write('lsbw,29'+'\n')
      f.write('wpoff,,,'+str(b_cc-a)+'\n')
      f.write('lsbw,24'+'\n')
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,,90'+'\n')
      f.write('wprota,,'+str((math.asin(B_R/C_R) + (a/C_R))*180/math.pi)+'\n')
      f.write('lsbw,10'+'\n')
      f.write('wprota,,'+str(((b_cs-a)/C_R)*180/math.pi)+'\n')   
      f.write('lsbw,24'+'\n')   
      f.write('wprota,,'+str(((b_cs-a)/C_R)*180/math.pi)+'\n') 
      f.write('lsbw,10'+'\n')
      f.write('WPCSYS,-1,0 \n') 
      f.write('wpoff,'+str(C_R)+',,'+str(B_R+(z))+'\n') 
      f.write('wpro,,-'+str(theta_cut)+'\n')
      f.write('asbw,9'+'\n') 
      f.write('WPCSYS,-1,0 \n')
      f.write('allsel,all'+'\n')
      f.write('wpoff,,,'+str(min(B_R+z+z_1,0.95*C_M_L))+'\n') 
      f.write('asbw,3'+'\n')
      f.write('asbw,8'+'\n') 
      f.write('WPCSYS,-1,0 \n')
      f.write('wpoff,,,'+str(B_R+z)+'\n') 
      f.write('asbw,1'+'\n')
      f.write('asel,u,,,6'+'\n') 
      f.write('asel,u,,,1'+'\n') 
      f.write('aatt,1,,1,0,1'+'\n')
      f.write('asel,s,,,6'+'\n') 
      f.write('asel,a,,,1'+'\n')
      f.write('aatt,1,,1,0,2'+'\n')
      f.write('asel,all'+'\n')    
      f.write('lsel,s,,,18'+'\n')
      f.write('lsel,a,,,32'+'\n')
      f.write('lesize,all,,,'+str(num_weld_c)+',,,,,0'+'\n')
      f.write('lsel,s,,,34'+'\n')
      f.write('lsel,a,,,30'+'\n')
      f.write('lsel,a,,,31'+'\n')
      f.write('lsel,a,,,33'+'\n')
      f.write('lesize,all,,,'+str(c)+',,,,,0'+'\n')
      f.write('lsel,s,,,29'+'\n')
      f.write('lsel,a,,,24'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut-(num_weld_c+2*c))+',,,,,0'+'\n')
      f.write('lsel,s,,,35'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut)+',,,,,0'+'\n')
      f.write('lsel,s,,,10'+'\n')
      f.write('lsel,a,,,19'+'\n')
      f.write('lsel,a,,,38'+'\n')
      f.write('lsel,a,,,4'+'\n')
      f.write('lesize,all,,,'+str(num_cir_crown)+',,,,,0'+'\n')
      f.write('lsel,s,,,36'+'\n')
      f.write('lsel,a,,,20'+'\n')
      f.write('lesize,all,,,'+str(num_cir_saddle)+',,,,,0'+'\n')
      f.write('allsel,all'+'\n')  
      f.write('lsel,s,,,27'+'\n')
      f.write('lsel,a,,,28'+'\n')
      f.write('lsel,a,,,17'+'\n')
      f.write('lsel,a,,,25'+'\n')
      f.write('lesize,all,,,'+str(d)+',,,,,0'+'\n')
      f.write('lsel,s,,,8'+'\n')
      f.write('lsel,a,,,16'+'\n')
      f.write('lesize,all,,,'+str(num_weld_b)+',,,,,0'+'\n')
      f.write('lsel,s,,,14'+'\n')
      f.write('lsel,a,,,26'+'\n')
      f.write('lesize,all,,,'+ str(num_br_cut-(num_weld_b+2*d))+',,,,,0'+'\n')
      f.write('allsel,all'+'\n')
      f.write('lsel,s,,,1'+'\n')
      f.write('lsel,a,,,2'+'\n')
      f.write('lesize,all,,,'+ str(e)+',,,,,0'+'\n')
      f.write('allsel,all'+'\n')  
      f.write('lsel,s,,,39'+'\n')
      f.write('lsel,a,,,40'+'\n')
      f.write('lsel,a,,,41'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut_1)+',,,,,0'+'\n')   
      f.write('lsel,s,,,13'+'\n')
      f.write('lsel,a,,,37'+'\n')
      f.write('lsel,a,,,11'+'\n')
      f.write('lesize,all,,,'+ str(num_chord)+',,,,,0'+'\n')
      f.write('lsel,s,,,6'+'\n')
      f.write('lsel,a,,,22'+'\n')
      f.write('lsel,a,,,12'+'\n')
      f.write('lsel,a,,,9'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_bt)+',,,,,0'+'\n')
      f.write('lsel,s,,,3'+'\n')
      f.write('lsel,a,,,5'+'\n')
      f.write('lesize,all,,,'+ str(num_brace)+',,,,,0'+'\n') 
      f.write('allsel,all'+'\n')
      f.write('asel,u,,,1'+'\n')
      f.write('asel,u,,,6'+'\n') 
      f.write('areverse,all'+'\n')  
      # Meshing 
      f.write('allsel,all'+'\n')
      f.write('AMAP,2,23,24,38,16'+'\n') 
      f.write('AMAP,5,17,38,24,4'+'\n') 
      f.write('AMAP,6,16,17,22,21'+'\n')
      f.write('amesh,4'+'\n')
      f.write('AMAP,10,23,39,40,24'+'\n') 
      f.write('AMAP,9,18,19,40,39'+'\n') 
      f.write('ADELE,7,,,1'+'\n') 
      f.write('ADELE,11,,,1'+'\n') 
      f.write('ADELE,3,,,1'+'\n')  
      f.write('AMAP,1,21,22,15,14'+'\n') 
      #to get node numbering file 
      f.write('ksel,s,,,26,37,1'+'\n')
      f.write('nslk'+'\n')  
      f.write('''nwrite,nodeindex,'txt',,'''+'\n')
      #Axial load
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('DSYM,SYMM,X, , '+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('DSYM,SYMM,Z,,'+'\n') 
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('DSYM,SYMM,Y,,'+'\n')     
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,UX,UY,UZ , , ,'+'\n')
      f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
      f.write("SF,all,pres,-"+ str(B_T)+'\n')
      #Inplane Bending
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('DSYM,SYMM,X, , '+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('DSYM,ASYM,Z,,'+'\n') 
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('DSYM,SYMM,Y,,'+'\n')     
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,ALL,, , , '+'\n') 
      div = 4
      for i in range(div):
          f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
          f.write("nsel,r,loc,z,"+ str(i*B_R/div)+","+str( 0.01+(i+1)*B_R/div)+'\n')
          f.write("SF,all,pres,-"+ str(B_T*(i+0.5)/div)+'\n')
      #Outplane Bending
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('DSYM,ASYM,X, , '+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('DSYM,SYMM,Z,,'+'\n') 
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('DSYM,SYMM,Y,,'+'\n')     
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,ALL,, , , '+'\n') 
      div = 4
      for i in range(div):
          f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
          f.write("nsel,r,loc,x,"+ str(i*B_R/div)+","+str( 0.01+(i+1)*B_R/div)+'\n')
          f.write("SF,all,pres,-"+ str(B_T*(i+0.5)/div)+'\n')
      #analysis    
      f.write('allsel,all'+'\n')  
      f.write('/sol'+'\n')
      f.write('antype,0'+'\n')
      f.write('solve'+'\n') 
      f.write('save\n')
      f.write('finish\n')
      if count<len(df)-1:
            f.write('BOPTN,NWARN,1\n')
            f.write('/UIS,MSGPOP,3 \n') 
            f.write('KEYW,PR_SGVOF,1  \n')
            f.write('/NERR,5,10000, ,1,5, \n')
            f.write('/clear \n')
            f.write('''/CWD,'''+"'"+directory+ str(count+1)+"'"+ '\n')
            f.write('''/filnam, '''+ str(count+1)+'\n')
            f.write('''/title, '''+ str(count+1)+'\n')
            f.write('''/INPUT,'model','txt','''+"'"+directory+ str(count+1)+'''\\',, 0'''+'\n')   
# saving the dataframe
str_folder = str(folder)+'_without_result'+'.csv'
df.to_csv(str_folder)

