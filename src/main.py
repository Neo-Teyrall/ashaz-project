import copy
import datetime

def creat_mat(seq):
    out  = []
    out.append(copy.deepcopy(seq))
    for i in range(len(seq)-1):
        seq = [seq[-1]] + seq[:-1]
        out.append(copy.deepcopy(seq))
    return out


def get_trasforme(seq):
    print("start transform")

    seq = "$" +  seq
    seq = list(seq)
    mat = creat_mat(seq)
    mat.sort()
    T = []
    c = 0
    cp = len(mat)
    for i in mat:
        T.append(i[-1])
        print(c,"/",cp,end = "\r")
        c += 1 
    return T

def get_chain(T) -> list:
    index = []
    c_i = {"A" : 0,
           "T" : 0,
           "C" : 0,
           "G" : 0,
           "$" : 0}

    first_pos = {"A":-1, "T" : -1, "C": -1, "G": -1 ,"$" : -1}
    
    sT = copy.deepcopy(T)
    sT.sort()

    for i in range(len(T)):
        for k in first_pos.keys():
            if sT[i] == k:
                if  first_pos[k] == -1:
                    first_pos[k] = i
            if T[i] == k : 
                index.append(c_i[k])
                c_i[k] += 1
                
    f_mat = [sT, T , index, index]
   
    pos = 0
    out = ""

    poss = []
    for i in range(len(T)):
        out = f_mat[1][pos] + out
        pos2 =first_pos[f_mat[1][pos]] + f_mat[2][pos] 
        f_mat[3][pos] = len(T)-i-1
        pos = pos2
        pass
    out = out[1:]
    return(out,f_mat[3])

def get_pos(seq):
    T = get_trasforme(seq)
    sort_T = copy.deepcopy(T)
    sort_T.sort()
    index = []
    mat_counter_i, index = get_mat_counter(T)
    first_pos = get_first_pos(T)


    return(get_position(T,first_pos,index))


def creat_counter(T,val = 0) :
    out = {}
    for i in T:
        if i not in out :
            out.setdefault(i,val)
    return out
        
def get_first_pos(T):
    sT = copy.deepcopy(T)
    sT.sort()
    out = creat_counter(T,-1)
    for i in range(len(T)):
        for k in out.keys():
            if sT[i] == k:
                if  out[k] == -1:
                    out[k] = i
    return out
def get_mat_counter(T):
    out = []
    index = []
    counter_i = creat_counter(T)
    for i in range(len(T)):
        for k in counter_i.keys():
            if T[i] == k :
                index.append(counter_i[k])
                counter_i[k] += 1
        out.append(copy.deepcopy(counter_i))
    return(out,index)

def find_read(T,sort_T,read,mat_counter_i,index,first_pos,pos):
    read = read[::-1]
    b , e = get_intervale_read(read, mat_counter_i, first_pos)
    if b > e :
        return None
    elif b == e :
        return(pos[b]+1)
    else :
        return([pos[b]+1,pos[e]+1])

def get_position(T,first_pos,index):
    position = [0]*len(T)
    pos = 0
    for i in range(len(T)):
        position[pos] = len(T)-i-2
        pos2 = first_pos[T[pos]] + index[pos] 
        pos = pos2
    return position

def get_intervale_read(read,mat_counter_i,first_pos):
    b : int = first_pos[read[0]]
    e : int = first_pos[read[0]] + mat_counter_i[-1][read[0]] -1
    for i in read[1:]:
        b = first_pos[i] + mat_counter_i[b-1][i]
        e = first_pos[i] + mat_counter_i[e][i] - 1
    return(b,e)
    pass

import sys
from pytictoc import TicToc
import char
from char import i

###############################################################################
################################## Code projet ################################
###############################################################################


def get_T(seq) -> list:
    index = [] # creation d'une liste
    # remplissage de la liste avec une liste [posiiton , character]
    for i, char in enumerate(seq):
        if i == 0 : # pour la premiere position prendre le dernier charracter
            c = seq[-1]
        else :
            c = seq[i-1] # pour les autre le caratere a la position antérieur
        add = [i,c]
        index.append(add) # ajouté a la liste
        
    for i, char_i in enumerate(seq):
        print("{:.2f} %".format(i/len(seq)*100),end = "\r") # affichage de progression
        for j, char_j in enumerate(seq[i+1:]):
            i2 = i
            j2 = j+i+1
            while j2+1 < len(seq)  and seq[i2] != seq[j2]  and i2 < j:
                i2 += 1
                j2 += 1
            if seq[i2] < seq[j2] and index[i2][0] > index[j2][0]:
                index[i2][0], index[j2][0] = index[j2][0], index[i2][0]
            elif seq[i2] > seq[j2] and index[i2][0] < index[j2][0]:
                index[j2][0], index[i2][0] = index[i2][0], index[j2][0]
            elif seq[i2] == seq[j2] and index[i2][0] < index[j2][0] and j2 == len(seq)-1:
                index[j2][0], index[i2][0] = index[i2][0], index[j2][0]


    index.sort()
    # récupéré la transformer.
    for i,val in enumerate(index) :
        index[i] = val[1]
    return index



if __name__ == "__main__" :
    in_seq= {}
    ok = False
    key = ""
    with open("./Hu-1.fasta.bwt") as filin:
        lines = filin.readlines()
        T = list(lines[0][:-1])
        pos = list(lines[1][:-2].split(" "))
        for i, p  in enumerate(pos):
            pos[i] = int(p)

    ok = False
    key = ""
    in_read = {}
    with open("./READSsars_cov_2_1e6.fasta") as filin:
        lines = filin.readlines()
        for i ,line in enumerate(lines):
            if line[0] == ">" :
                if not ok :
                    ok = True
                else:
                    in_read.setdefault(key,string)
                key = line[:-1]
                string = ""
            else :
                string += line[:-1]
        in_read.setdefault(key,string)

    out1 = {}
    out2 = {}
    c = 0
    ct = len(in_read)
    get_T()
    from time import time
    t1 = time()
    # do something
    t2 = time()
    sort_T = copy.deepcopy(T)
    sort_T.sort()
    mat_counter_i, index = get_mat_counter(T)
    first_pos = get_first_pos(T)
    pos = get_position(T,first_pos,index)
    for j in in_read:
        s = (t2 - t1) * (ct-c)
        m = s/ 60
        H = int(m/60)
        #print("{:2}:{:2}:{:2}".format(H,int(m%60),int(s%60)), "\t seconde",end="\r")
        print("{:.2f}".format(c/ct*100),"%",end = "\r")
        t1 = time()
        c +=1
        t.tic()
        a = find_read(T,sort_T, in_read[j],mat_counter_i,index,first_pos,pos)
        if a== None:
            t2 = time()
            continue
        elif isinstance(a,list):
            out2.setdefault(j,a)
        else:
            out1.setdefault(j,a)
        t2 = time()
    print(out1)
    print(out2)
    seq = ""
    with open("./Hu-1.fasta") as filin:
        lines = filin.readlines()
        for i, line in enumerate(lines):
            if line[0] == ">":
                continue
            seq += line[:-1] 
    a = get_T(seq)

            
            
    
