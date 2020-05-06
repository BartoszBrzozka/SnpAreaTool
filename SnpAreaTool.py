from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from Bio import SeqIO
import gzip
import numpy as np
import os

 
class Root(Tk):
    def __init__(self):
        super(Root, self).__init__()
        self.title("SNP Area Tool")
        self.minsize(300, 310)

        self.labelFrame = ttk.LabelFrame(self, text = "Open VCF/Bed file")
        self.labelFrame.grid(column = 0, row = 4, padx = 1, columnspan = 3, pady = 1, sticky = W)
 
        self.button1()

        self.labelFrame1 = ttk.LabelFrame(self, text = "Open reference - Fasta file")
        self.labelFrame1.grid(column = 0, row = 5, columnspan = 3, padx = 1, pady = 1,sticky = W)
 
        self.button2()

        ttk.Label(self, text="<-- bp").grid(row=0,column = 1, pady=1)
        ttk.Label(self, text="bp -->").grid(row=1,column =1, pady=1)
        global d1,d2
        d1 = StringVar()
        d2 = StringVar()
        e1 = ttk.Entry(self,textvariable=d1)
        e2 = ttk.Entry(self,textvariable=d2)
        e1.grid(row=0, column=2, pady=1, sticky = W)
        e2.grid(row=1, column=2, pady=1, sticky = W)
        
        self.labelFrame2 = ttk.LabelFrame(self)
        self.labelFrame2.grid(column = 4, row = 0,rowspan =2, padx = 10, pady = 0, sticky = E)

        self.listbox()
        self.listbox.grid(row=6,column=0,padx=5,pady=1,columnspan =5)

    def button1(self):
        self.button = ttk.Button(self.labelFrame, text = "Wybierz plik", command = self.fileDialog)
        self.button.grid(column = 1, row = 1)
        
    def button2(self):
        self.button = ttk.Button(self.labelFrame1, text = "Wybierz plik", command = self.fileDialog2)
        self.button.grid(column = 1, row = 1)
        
    def listbox(self):
        self.listbox= Listbox(self,width = 50)

## funkcja pobierająca koordynaty z entryboxów ##
    def getentry(self):
        global left,right
        self.listbox.delete('0','end')
        left = d1.get()
        right = d2.get()
        self.listbox.insert(END, "Coordinates " + left + " : " + right + " loaded")



## Funkcja pobierająca dane z plików VCF/BED ##
    def fileDialog(self):
        global plik
        plik = '' 
        self.listbox.delete('0','end') 
        if not ((d1.get().isnumeric() and d2.get().isnumeric())):
            self.listbox.insert(END,"Enter valid coordinates (for BED file enter 0)")
        else:
            if not (int(d1.get() and d2.get())<0):
                self.getentry()
                plik = filedialog.askopenfilename(initialdir =  "/", title = "Wybierz plik", filetype =
                (("VCF","*.vcf"),("BED","*.bed")))
                if not(plik == ''):
                    self.listbox.insert(END,os.path.basename(plik) + " file loaded")
                    if plik.endswith('.vcf'):
                        with open(plik, "r") as datafile:
                            for line in datafile:
                                li = line.strip()
                                if not (li.startswith("#")):
                                    fields = line.strip().split()
                                    if (fields[0].startswith("Chr")):
                                        chrom.append(fields[0].replace("Chr",""))
                                        positions1.append(int(fields[1])-int(left))
                                        positions2.append(int(fields[1])+int(right))
                    if plik.endswith('.bed'):
                        with open(plik, "r") as datafile:
                            for line in datafile:
                                fields = line.strip().split()
                                chrom.append(fields[0])
                                positions1.append(int(fields[1]))
                                positions2.append(int(fields[2]))    

                else:
                    self.listbox.insert(END,"File loading error")
               
## Funkcja pobierająca dane z plików fasta (możliwe pobranie danych ze spakowanych plików) ##
    def fileDialog2(self):
        global handle
        if (plik == ''):
            self.listbox.insert(END,"Choose VCF/BED file")
        else:
            self.filename = filedialog.askopenfilename(initialdir =  "/", title = "Wybierz plik", filetype =
            (("fa","*.fa"),("fasta","*.fasta"),("fa.gz","*.fa.gz"),("fasta.gz","*.fasta.gz")))

            self.listbox.insert(END,os.path.basename(self.filename) + " file loaded. Calculations starts...")
            
            if self.filename.endswith('.fa.gz'):
                with gzip.open(self.filename, "rt") as handle:
                    self.listbox.insert(END,self.filename + " file loaded")
                    self.calculate()
            elif self.filename.endswith('.fasta.gz'):
                with gzip.open(self.filename, "rt") as handle:
                    self.listbox.insert(END,self.filename + " file loaded")
                    self.calculate()
                        
            elif self.filename.endswith('.fa'):
                with open(self.filename, "r") as handle:
                    self.listbox.insert(END,self.filename + " file loaded")
                    self.calculate()

            elif self.filename.endswith('.fasta'):
                with open(self.filename, "r") as handle:
                    self.listbox.insert(END,self.filename + " file loaded")
                    self.calculate()

## funkcja porównująca dane z obu plików wejściowych ## 
    def calculate(self):
        print(len(chrom),len(positions1),len(positions2))
        for record in SeqIO.parse(handle, "fasta"):
            if "GK" in record.name:
                for j,k,l in zip(chrom,positions1,positions2):
                    if(record.description.split()[2] == j):
                        tablica.append([j,k,l,str(record.seq[k:l])])


        self.listbox.insert(END,"Sequences loaded")
        with open('tablica.txt', 'w') as f:
            for item in tablica:
               f.write("%s\n" % item)
        if plik.endswith('.vcf'):
            with open('tablica.txt', 'r') as f:
                wynik.extend(["CHROM","POS","POS-","POS+","A","G","C","T","N","%GC","%N"])
                wynik.extend("\n")
                lines = f.readlines()
                for k in lines:
                    line = k.strip().split()
                    c=g=a=t=n=0
                    for i in str(line[3]) :
                        if(i == "C"): c+=1
                        elif(i == "G"): g+=1
                        elif(i == "T"): t+=1
                        elif(i == "A"): a+=1
                        elif(i == "N"): n+=1
                    
                    if((a+g+c+t+n) != 0):
                        wynik.extend([line[0].replace("[","").replace("'","").replace(',',''),int(line[2].replace(',',''))-100,line[1].replace(',',''),line[2].replace(',',''),a,g,c,t,n,"{:.2f}".format(((c+g)/(a+g+c+t+n))*100),"{:.2f}".format(100*(n/(int(line[2].strip(','))-int(line[1].strip(',')))))])
                        wynik.extend("\n")
                    else:
                        wynik.extend([line[0].replace("[","").replace("'","").replace(',',''),int(line[2].replace(',',''))-100,line[1].replace(',',''),line[2].replace(',',''),a,g,c,t,n,0.0,"{:.2f}".format(100*(n/(int(line[2].strip(',')-int(line[1].strip(','))))))])
                        wynik.extend("\n")
                
                self.filename = filedialog.asksaveasfilename(initialdir =  "/", title = "Save file")
                with open(self.filename, 'w') as f:
                    for item in wynik:
                        if(item != "\n"):
                            f.write("%s\t" % item)
                        else:
                            f.write("\n")
                self.listbox.insert(END,"File created")
                wynik.clear()
                tablica.clear()
                chrom.clear()
                positions1.clear()
                positions2.clear()
        elif plik.endswith('.bed'):
            with open('tablica.txt', 'r') as f:
                wynik.extend(["CHROM","POS-","POS+","A","G","C","T","N","%GC","%N"])
                wynik.extend("\n")
                lines = f.readlines()
                for k in lines:
                    line = k.strip().split()
                    c=g=a=t=n=0
                    for i in str(line[3]) :
                        if(i == "C"): c+=1
                        elif(i == "G"): g+=1
                        elif(i == "T"): t+=1
                        elif(i == "A"): a+=1
                        elif(i == "N"): n+=1
                    
                    if((a+g+c+t+n) != 0):
                        wynik.extend([line[0].replace("[","").replace("'","").replace(',',''),line[1].replace(',',''),line[2].replace(',',''),a,g,c,t,n,"{:.2f}".format(((c+g)/(a+g+c+t+n))*100),"{:.2f}".format(100*(n/(int(line[2].replace(',',''))-int(line[1].replace(',','')))))])
                        wynik.extend("\n")
                    else:
                        wynik.extend([line[0].replace("[","").replace("'","").replace(',',''),line[1].replace(',',''),line[2].replace(',',''),a,g,c,t,n,0.0,"{:.2f}".format((n/(int(line[2].replace(',',''))-int(line[1].replace(',','')))))*100])
                        wynik.extend("\n")
        
            self.filename = filedialog.asksaveasfilename(initialdir =  "/", title = "Save file")
            with open(self.filename, 'w') as f:
                for item in wynik:
                    if(item != "\n"):
                        f.write("%s\t" % item)
                    else:
                        f.write("\n")
            self.listbox.insert(END,"File created")
            wynik.clear()
            tablica.clear()
            chrom.clear()
            positions1.clear()
            positions2.clear()
            


global left
global right
global handle  
global plik
wynik = []
tablica = []
chrom = []
positions1 = []
positions2 = []
root = Root()
root.mainloop()
