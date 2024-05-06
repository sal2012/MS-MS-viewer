import xml.etree.ElementTree as ET
import gzip
import sys
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt



filename = sys.argv[1]
scan_number = sys.argv[2]
peptide_sequence = sys.argv[3]
 

ns = '{http://sashimi.sourceforge.net/schema/}'

# I made this for loop to parse xml file where interparse read xml file element by element.
for evt,ele in ET.iterparse(filename):
    
    #this if statement will chech if this element have this tag scan
	if ele.tag == ns + 'scan':
        #print(ele.attrib) #this prints all of the scan w the num
            

        #if it finds the word scan+ns
        #then we r gonna check if scan tag have attribute called number that is=scan_num            
		if ele.attrib.get('num') == scan_number:
	    
            #print(ele.attrib)
                    
#if i find the scan no i want from the if the statement,
#the following for loop isolates the peaks element
                    #where it look at every element to find the word peaks
			for elem in ele:
                            
        #if it find the peaks inside the element    
				if elem.tag == ns + 'peaks':
                                    
                    #print(elem)

                    #saving the peaks data into a variable called peakselt                  
					peakselt = elem
                    #print(peakselt)
                    
        
# peakselt is the XML element corresponding to the peaks list
peaks = array('f',b64decode(peakselt.text))
if sys.byteorder != 'big':
    peaks.byteswap()
mzs = peaks[::2]
ints = peaks[1::2]

#print(mzs)
#print(ints)


#secondstep
# I made a dictionary for the masses of the amino acids
 
aa_monoisotopic_mass = {
    'A': 71.04, 
    'C': 103.01, 
    'D': 115.03, 
    'E': 129.04, 
    'F': 147.07,
    'G': 57.02, 
    'H': 137.06, 
    'I': 113.08, 
    'K': 128.09, 
    'L': 113.08,
    'M': 131.04, 
    'N': 114.04, 
    'P': 97.05, 
    'Q': 128.06, 
    'R': 156.10,
    'S': 87.03, 
    'T': 101.05, 
    'V': 99.07, 
    'W': 186.08, 
    'Y': 163.06
}

   

#I initialized these 2 dictionaries 2 store all the values of b ions and y ions that i am gonna get 
b_ions = {}
y_ions = {}
#this the proton mass
b_mass = 1
#this the water mass
y_mass = 19  

# I Calculated the b-ions by iterationg through peptide sequence to calculate the m/z values of b ions in the sequence except the last one I will add its mass to b_mass and store mass plus the proton mass in the b ions list.
count = 0

#I used range bec it will produce seq of numbers that is every pos in peptide seq
#and I put -1 bec I don't won't whole lenght
for i in range(len(peptide_sequence) - 1):
    count += 1

    #it will access dic for aa masses
    #adds and assigns the mass of AA at pos i in peptide_sequence to b-mass .
    b_mass += aa_monoisotopic_mass[peptide_sequence[i]]
    #so I made a dict called b_ions where key is B with count and the value is the mass I calc
    b_ions['B' + str(count)] = (b_mass)

# ICalculated y-ion
count = 0

#-1 start from the back of the peptide seq
    #, and move backwards
for i in range(len(peptide_sequence) - 1, 0, -1):
    count += 1
        
    #it will access dic for aa masses
    #adds the mass of AA at pos i in peptide_sequence to y_mass.
    y_mass += aa_monoisotopic_mass[peptide_sequence[i]]

    #i made akey for dic y_ions,It combines y w current count converted to string.
    #if count is 1, the key will be y1= y_mass(value) in the b_ions dictionary.
    y_ions['Y' + str(count)] = (y_mass)  


# i created alist that have these 2 elements and put it invariable mz_int
mz_int = list(zip(mzs,ints))
tolerance = 0.2
matched_b_ions = []
matched_y_ions = []
#I Matched observed peaks to calculate b-ion and y-ion m/z values

# Matching for b-ions

#i made for loop to acces mz and intensity in each iteration
for mz, intensity in mz_int:
    
    # ia made another for loop inside it to iterate in b ion dict where ion_label is the key
    #and b_mz is the value
    for ion_label, b_mz in b_ions.items():
        
        #i did the if statement bec I want to check if observed mz is within
        # 0.2 toleranceof the calculated m/z value 
        if abs(mz - b_mz) <= tolerance:
            
            # I made this list to put the matched ion and it's m/z and intensty
            matched_b_ions.append((mz, intensity, ion_label))

# Matching for y-ions
#i made for loop to acces mz and intensity in each iteration
for mz, intensity in mz_int:

    # ia made another for loop inside it to iterate in y ion dict where ion_label is the key
    #and y_mz is the value
    for ion_label, y_mz in y_ions.items():

        #i did the if statement bec I want to check if observed mz is within
        # 0.2 toleranceof the calculated m/z value         
        if abs(mz - y_mz) <= tolerance:

        # I made this list to put the matched ion and it's m/z and intensty    
            matched_y_ions.append((mz, intensity, ion_label))


#plotting

#i set x is mzs and y intensity, the stem line color
plt.stem(mzs, ints, linefmt='grey',markerfmt=' ')
#plt.stem(mzs, ints, linefmt='grey')

#I made a for loop it extreact mz,intensity, ion label from matched_b_ ions
for mz, intensity, ion_label in matched_b_ions:

    #this creates stem plot for every matched b ion 
    plt.stem([mz], [intensity], linefmt='green', markerfmt=' ')

#thiswill add the text annotation at location of each matched b ions
# the mz, intensity to show the coordination on plot where it will be placed
for mz, intensity, ion_label in matched_b_ions:

      #this creates stem plot for every matched y ion 
    plt.text(mz, intensity, ion_label, color='green')

#I made a for loop it extreact mz,intensity, ion label from matched_ y_ions
for mz, intensity, ion_label in matched_y_ions:
    plt.stem([mz], [intensity], linefmt='purple', markerfmt=' ')

#thiswill add the text annotation at location of each matched y ions
# the mz, intensity to show the coordination on plot where it will be placed
for mz, intensity, ion_label in matched_y_ions:

    #this creates stem plot for every matched y ion    
    plt.text(mz, intensity, ion_label, color='purple')


plt.title("Mass Spectrometry Data - Matched b-ions and y-ions")
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.show()

