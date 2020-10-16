import random
dic = ["A","T","C","G"]
def gen(_id,x,f):
    out =""
    for i in range(x):
         out += dic[random.randint(0,3)]
    
    f.write("@"+str(_id)+"SEQ_ID\n")
    f.write(out+"\n")

def gen_set(nb_seq, taille, name):
    with open(name + ".fastq","w") as filout:
        for i in range(nb_seq):
            gen( i,taille, filout)

if __name__ == "__main__" :
    gen_set(100,5,"read")
    gen_set(100,100,"seq")
    pass
