import copy
import datetime
import sys
import argparse

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', metavar='FASTA FILE', nargs='?', help='fichier fasta du génome analysé')
    parser.add_argument('-r', metavar='READ FASTQ FILE', nargs='?', help='fichier fastq des séquences read recherchées')
    parser.add_argument('-o',   metavar='CSV FILE', help='fichier d\'output', default="results.csv")
    parser.add_argument('-m', default = 1, type=int, help="Recherche read sans mutation, mettre à 0 pour désactiver l'option.")
    parser.add_argument('-i', default = 0, type=int, help="Recherche read avec insertion, mettre à 1 pour activer l'option.")
    parser.add_argument('-d', default = 0, type=int, help="Recherche read avec deletion, mettre à 1 pour activer l'option.")
    parser.add_argument('-s', default = 0, type=int, help="Recherche read avec substitution, mettre à 1 pour activer l'option.")
    parser.add_argument('-l', default = 1, type=int, help="Taille des insertion/deletion/substitution, mettre à 1 pour activer l'option.")
    args = parser.parse_args()

    if args.f and args.f[-5:] != "fasta":
        raise ValueError("File {} is not FASTA".format(args.f))

    if args.r and args.r[-5:] != "fastq":
        raise ValueError("File {} is not FASTQ".format(args.r))

    return args

def read_fastq(fastq):
    """Lit le fichier <fastq>.fastq passé en paramètre.

    Parametres
    ---------
    fastq : fichier fastq
    ------
    reference: str
        la séquence de référence
    """
    with open(fastq, "r") as filin:
        for line in filin:
            yield next(filin)[:-1]


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


def creat_counter(T,val = 0) :
    out = {}
    for i in T:
        if i not in out :
            out.setdefault(i,val)
    return out

def find_read(read,mat_counter_i,first_pos,pos):
    """
    Si b > e alors le read ne correspond à aucune position, s'ils sont égaux le read est à la position indiquée par eux. 
    Si e > b il y a plusieurs poisitions au read dont b et e.
    """
    b , e = get_intervale_read(read, mat_counter_i, first_pos)
    if b > e :
        return None
    elif b == e :
        return([pos[b]+1])
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



################################################################################
####################                get read                ####################
################################################################################


def get_intervale_read(read,mat_counter_i,first_pos):
    """
    Fonction prenant en argument d'entrée le read, le compteur et la première position.
    Elle permet de retourner b et e, les positions potentielles du read sur la seq.
    """
    read = read[::-1]
    b : int = first_pos[read[0]]
    e : int = first_pos[read[0]] + mat_counter_i[-1][read[0]] -1
    for i in read[1:]:
        b = first_pos[i] + mat_counter_i[b-1][i]
        e = first_pos[i] + mat_counter_i[e][i] - 1
    return(b,e)


def get_intervale_read_sub(in_read, mat_counter_i,
                           first_pos, subs = 1):
    """
    Fonction prenant en argument d'entrée le read, le compteur, la première position et le nombre de substitutions acceptées.
    Elle permet de retourner b et e, les positions potentielles du read sur la seq.
    """
    # initialisation des paraetre
    read = in_read[::-1] 
    AN = ["A","T","C","G"]
    pe = -1
    pb = -1
    sub_pos = -1
    subs_an = []
    i = 1
    nb_sub = 0
    ## debut
    b : int = first_pos[read[0]]
    e : int = first_pos[read[0]] + mat_counter_i[-1][read[0]] -1
    
    while i < len(read):
        an = read[i] # l'acide aminé actuelle
        pb , pe = b, e # sauvegardé les anciennes position.
        b = first_pos[an] + mat_counter_i[b-1][an]
        e = first_pos[an] + mat_counter_i[e][an] - 1
        # si le read continure 
        if i< len(read) and b <= e:
            i += 1
            print(i-1 ,"->", i )
        # erreure dans le read
        if b > e:
            print(read,"to sub",read[i],i)
            if sub_pos != i : # la position n'a pas encore été substituer 
                sub_pos = i
                subs_an = []
                nb_sub += 1
                # si le nombre de substitution est trop grand
                if nb_sub > subs:
                    print("trop de substitution")
                    return(b,e)
            # substituer la position #
            #// ajouter l'acide nucleique dans la liste des acide nucléique déja
            # tester \\#
            subs_an.append(read[i])
            # choisir un acide aminé dans les acide nucélique  pas encore tester
            new_an = list(set(AN)-set(subs_an))[0]
            # remplacer l'acide nucléique
            read = read[:i] + new_an + read[i+1:]
            ## retour au valeur de b et e de la position antétieru 
            e = pe
            b = pb
            print(read,"subted",read[i])
    return(b,e)

def get_intervale_read_del(in_read, mat_counter_i, first_pos,
                           dels = 1):
    """
    Fonction prenant en argument d'entrée le read, le compteur, la première position et le nombre de délétions acceptées.
    Elle permet de retourner b et e, les positions potentielles du read sur la seq.
    """
    # initialisation des variable.
    in_read_reverse = in_read[::-1]
    read =in_read_reverse
    pe = -1
    pb = -1
    sub_pos = -1
    subs_an = []
    i = 1
    nb_dels = 0
    b : int = first_pos[read[0]]
    e : int = first_pos[read[0]] + mat_counter_i[-1][read[0]] -1
    
    while i < len(read):
        an = read[i]
        pb , pe = b, e
        b = first_pos[an] + mat_counter_i[b-1][an]
        e = first_pos[an] + mat_counter_i[e][an] - 1
        # si la il ny a pas d'erreure dans l'alignement 
        if i< len(read) and b <= e:
            i += 1
            print(i-1 ,"->", i )
        # si l'aligenement est incorecte.
        if b > e:
            print(read,"to del",read[i],i)
            nb_dels += 1 # augementé le compteur de déléttion 
            
            if nb_dels > dels: # Si trops 
                print("trop de substitution")
                return(b,e)
            # suprimé la postion qui ne s'alligne pas
            read = read[:i] + read[i+1:]
            # reprendre la b et e de la position anterieurs.
            e = pe
            b = pb
            print(read,"delated")
    return(b,e)



def get_intervale_read_ins(in_read, mat_counter_i, first_pos,
                           inss = 1):
    """
    Fonction prenant en argument d'entrée le read, le compteur, la première position et le nombre d'insertions acceptées.
    Elle permet de retourner b et e, les positions potentielles du read sur la seq.
    """
    in_read_reverse = in_read[::-1]
    read =in_read_reverse
    AN = ["A","T","C","G"]
    pe = -1
    pb = -1
    ins_pos = -1
    ins_an = []
    b : int = first_pos[read[0]]
    e : int = first_pos[read[0]] + mat_counter_i[-1][read[0]] -1
    i = 1
    nb_ins = 0
    while i < len(read):
        an = read[i]
        pb , pe = b, e
        b = first_pos[an] + mat_counter_i[b-1][an]
        e = first_pos[an] + mat_counter_i[e][an] - 1
        if i< len(read) and b <= e:
            i += 1
            print(i-1 ,"->", i )
        if b > e:
            if ins_pos != i :
                print(read,"to to ins", i-1, i,)
                ins_pos = i
                ins_an = []
                nb_ins += 1

                if nb_ins > inss:
                    print("trop d'insertion")
                    return(b,e)
                new_an = list(set(AN)-set(ins_an))[0]
                ins_an.append(new_an)
                read = read[:i] + new_an + read[i:]
                e = pe
                b = pb
            else :
                print(read,"to sub", read[i], i)
                new_an = list(set(AN)-set(ins_an))[0]
                ins_an.append(new_an)
                read = read[:i] + new_an + read[i+1:]
                e = pe
                b = pb
            print(read,"insered",read[i])
    return(b,e)

###############################################################################
################################## Code projet ################################
###############################################################################

def get_seq(file_path : str):
    """
    Fonction prenant en argument le chemin vers le fichier fasta contenant la séquence à analyser et la récupère. 
    """
    seq = ""
    with open(file_path) as filin:
        lines = filin.readlines()
        for i, line in enumerate(lines):
            if line[0] == ">":
                continue
            seq += line[:-1]
    return seq


def get_T(seq:str) -> list:
    """
    Fonction prenant la sequénce en argument. Elle va créer la transformée de Burrows-Wheeler en liant chaque nucléotide de la séquence à une position. 
    Puis en inversant la position de 2 nucléotides si le nucléotide suivant ou la série suivante de nucléotide arrive plus tôt dans l'ordre alphabétique alors que la position est plus élevée.
    """
    index = [] # creation d'une liste
    seq +=  "$"
    # remplissage de la liste avec une liste [posiiton , character]
    for i, char in enumerate(seq):
        if i == 0 : # pour la premiere position prendre le dernier caractere
            c = seq[-1]
        else :
            c = seq[i-1] # pour les autre le caractere a la position antérieur
        add = [i,c]
        index.append(add) # ajouté a la liste
    index  = index[1:] + [index[0]]
    for i, char_i in enumerate(seq):
        print("{:.2f} %".format(i/len(seq)*100),end = "\r") # affichage de progression
        for j, char_j in enumerate(seq[i+1:]):
            i2 = i
            jb = i+j+1
            j2 = jb
            while j2+1 < len(seq)  and seq[i2] == seq[j2] :
                i2 += 1
                j2 += 1
            if seq[i2] < seq[j2] and index[i][0] > index[jb][0]:
                index[i][0], index[jb][0] = index[jb][0], index[i][0]
            elif seq[i2] > seq[j2] and index[i][0] < index[jb][0]:
                index[jb][0], index[i][0] = index[i][0], index[jb][0]
            elif seq[i2] == seq[j2] and index[i][0] < index[jb][0] and j2 == len(seq)-1:
                index[j2][0], index[i2][0] = index[i2][0], index[j2][0]

    for i in range(len(index)):
        if i  == 0 :
            last = index[-1][0]
            index[-1][0] = index[i][0]
        elif i== len(index)-1:
            index[i-1][0] = last
        else:
            index[i-1][0] = index[i][0]
    index.sort()
    for i,val in enumerate(index) :
        index[i] = val[1]
    return index



def get_mat_counter(T):
    """
    Prends en argument d'entrée la transformée et lui donne un compteur.
    """
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

def get_first_pos(T):
    """
    Donne la première position de la transformée
    """
    sT = copy.deepcopy(T)
    sT.sort()
    out = creat_counter(T,-1)
    for i in range(len(T)):
        for k in out.keys():
            if sT[i] == k:
                if  out[k] == -1:
                    out[k] = i
    return out


###############################################################################
##########                          Main                            ###########
###############################################################################


if __name__ == "__main__" :

    args = get_arguments()
    print(args)
    if(args.f):
        #Récupère la séquence
        seq = get_seq(args.f)

        #Calcul de la transformée
        T = get_T(seq)

        # Création du compteur 
        mat_counter_i,index = get_mat_counter(T)
    
        # Récupération de la première position
        first_pos  = get_first_pos(T)
    
        # Récupération de l'ensemble des positions
        pos = get_position(T,first_pos,index)

        with open(args.o, "w") as fillout:
            fillout.write("Read\tType\tPosition\tTaille arrangement\n")
            #Pour chaque read on va rechercher les matchs
            for read in read_fastq(args.r):
                if args.m == 1:
                    interval = find_read(read, mat_counter_i, first_pos, pos)
                    if interval != None:
                        for i in interval:
                            fillout.write("{}\tnormal\t{}\t0\n".format(read, i))

                #Insertion
                if args.i == 1:
                    b,e = get_intervale_read_ins(read, mat_counter_i, first_pos, inss = args.l)
                    if b==e:
                        fillout.write("{}\tInsertion\t{}\t{}\n".format(read, e, args.l))
                    if e > b:
                        fillout.write("{}\tInsertion\t{}\t{}\n".format(read, e, args.l))
                        fillout.write("{}\tInsertion\t{}\t{}\n".format(read, b, args.l))
                    #print(get_intervale_read_ins("AG",mat_counter_i,first_pos,inss = args.l))

                #Deletion
                if args.d == 1:
                    b,e = get_intervale_read_del(read, mat_counter_i, first_pos, dels = args.l)
                    if b==e:
                        fillout.write("{}\tDeletion\t{}\t{}\n".format(read, e, args.l))
                    if e > b:
                        fillout.write("{}\tDeletion\t{}\t{}\n".format(read, e, args.l))
                        fillout.write("{}\tDeletion\t{}\t{}\n".format(read, b, args.l))

                    #print(get_intervale_read_del("ACCG",mat_counter_i,first_pos,dels = args.l))

                #Substitution
                if args.s == 1:
                    b,e = get_intervale_read_sub(read, mat_counter_i, first_pos, subs = args.l)
                    if b==e:
                        fillout.write("{}\tSubstitution\t{}\t{}\n".format(read, e, args.l))
                    if e > b:
                        fillout.write("{}\tSubstitution\t{}\t{}\n".format(read, e, args.l))
                        fillout.write("{}\tSubstitution\t{}\t{}\n".format(read, b, args.l))
                    #print(get_intervale_read_sub("TACG",mat_counter_i,first_pos, subs = args.l))




"""
if __name__ == "__main__" :
    in_seq= {}
    ok = False
    key = ""
    
    # Ouvre le fichier et récupère directement la transformée
    with open("./Hu-1.fasta.bwt") as filin:
        lines = filin.readlines()
        T = list(lines[0][:-1])
        pos = list(lines[1][:-2].split(" "))
        for i, p  in enumerate(pos):
            pos[i] = int(p)
    ok = False
    key = ""
    in_read = {}
    # Ouvre un fichier avec tous les reads
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
    
    # Récupère la séquence du fichier fasta
    seq = get_seq("./Hu-1.fasta")
    # Création de la séquence de manière brute 
    seq = "ATATCGT"
    # Création de la transformée de Burrows-Wheeler
    T = get_T(seq)

    # Création du compteur 
    mat_counter_i,index = get_mat_counter(T)
    
    # Récupération de la première position
    first_pos  = get_first_pos(T)
    
    # Récupération de l'ensemble des positions
    pos = get_position(T,first_pos,index)
    
    # Cherche la ou les positions d'un read précis
    print("yololo")
    interval = find_read("ATC", mat_counter_i, first_pos, pos)
    print(interval)
    print("bonjour")

    ### Les fonctions suivantes renvoies b et e directement. Il faut donc savoir que : si b > e alors le read ne correspond à aucune position, s'ils sont égaux le read est à la position indiquée par eux. Si e > b il y a plusieurs poisitions au read dont b et e.                     
    # Cherche la ou les positions d'un read précis                         
    print(get_intervale_read_sub("TACG",mat_counter_i,first_pos))
    # Cherche la ou les positions d'un read avec subs substitutions
    print(get_intervale_read_sub("TACG",mat_counter_i,first_pos, subs=1))
    # Cherche la ou les positions d'un read avec dels déléttions   
    print(get_intervale_read_del("ACCG",mat_counter_i,first_pos,dels = 2))
    # Cherche la ou les positions d'un read avec inss insertions    
    print(get_intervale_read_ins("AG",mat_counter_i,first_pos,inss = 1))

"""