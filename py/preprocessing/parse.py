import sys,os

fin = open("DiffCrossSections.h","r")
fout = open("DiffCrossSections2.h","w")

fxsec = open("updatedgasdiffxsec.txt","r")
xsec_new_v = []
for line in fxsec:
    words = line.split(",")
    for word in words:
        xsec_new_v.append(float(word))

energy_v = []
xsec_v = []

for line in fin:
    words = line.split(",")
    wordsnew = words
    print len(words)
    #if (len(words) < 10):
        #print line
    if (len(words) == 22510):
        ctr = -1
        xsecctr = 0
        newwords = []
        for word in words:
            ctr += 1
            try:
                xsecval = float(word)
                if ( (xsecval < 1e-12) and (xsecval > 0)):
                    #print 'energy : ',words[ctr-1]
                    #print 'xsec   : ',xsecval 
                    energy_v.append( float(words[ctr-1]) )
                    xsec_v.append(xsecval)
                    newwords.append(str(xsec_new_v[xsecctr]))
                    #words[ctr] = 
                    xsecctr += 1
                    #print xsecval
                else:
                    newwords.append(word)
            
            except:
                newwords.append(word)
                continue
        newstr = ""
        for i,word in enumerate(newwords):
            newstr += word
            if (i < len(newwords)-1):
                newstr += ","
            if (i == len(newwords)-1):
                newstr += '\n'
        fout.write(newstr)
            
    else:
        print line
        fout.write(line)

fin.close()
fout.close()
fxsec.close()

        
#print xsec_v
#print len(xsec_v)

#print xsec_new_v
#print len(xsec_new_v)


#print energy_v
#print len(energy_v)
