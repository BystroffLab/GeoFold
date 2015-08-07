#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import cgi
import urlparse
import os

import cgitb
cgitb.enable

def parseit(query, keyword):
    #print("<br>keyword: "+keyword+"<br>query: "+query+"<br>")
    key_index = query.find(keyword)
    #print(key_index)
    if key_index != -1:
        #print("<br>KEYINDEX<br>")
        andindex = query.find("&",key_index)
        #print(andindex)
        output = query[key_index+len(keyword)+1:andindex]
        #print(output)
        return output
    return ''
    
def readDagFile(iseg,dag):
    iststate = False
    try:
        readDag = open(dag,'r')
    except IOError:
        oneliner("DAG file is missing: "+dag)
        readDag.close()
        return
    if iseg[0] == 'n': #intermediate
        try:
            inode = int(iseg[1:])
        except ValueError:
            oneliner("Bad iseg number: "+iseg)
            readDag.close()
            return
        iststate = False
    elif iseg[0] == 't': #transition state
        try:
            inode = int(iseg[2:])
        except ValueError:
            oneliner("Bad iseg number: "+iseg)
            readDag.close()
            return
        iststate = True
    else:
        oneliner("ERROR bad node label: "+iseg)
        return
    while 1:
        line = readDag.readline()
        if line == '':
            return
        if line[0:6] == "TSTATE" and iststate:
            try:
                i = int(line[6:14])
            except ValueError:
                oneliner("Error: Bad i number"+line[6:14])
                return
            if i == inode:
                readDag.close()
                tstateline(line,dag)
                return
        if line[0:6] == "ISEGMT" and not iststate:
            try:
                i = int(line[6:14])
            except ValueError:
                print line
                oneliner("Error: Bad i number"+line[6:14])
                return
            if i == inode:
                bline = readDag.readline()
                readDag.close()
                istateline(line,bline,dag)
                return
                
def tstateline(line,dag):
    cuttypes = {'h':"HINGE",'p':"PIVOT",'b':"BREAK",'m':"MELTING",'s':"SEAM",'u':"UNKNOWN CUTTYPE"}
    line = line.split()
    nts = line[1]
    f = line[2]
    u1 = line[3]
    u2 = line[4]
    ntrp = line[5]
    ctype = line[6]
    iseam = line[7]
    traf = line[8]
    cuttype = cuttypes[ctype]
    print '<font color="#000055"><pre>%s %s %s %s %s %s %s %s %s</pre></font>'%tuple(line)
    print '<br><b>Transition state %9i </b>&nbsp;%s'%(int(nts),cuttype)
    print '<br>folded state node number = <a href="isegment.cgi?iseg=n%s&dag=%s&">%s</a>'%(f,dag,f)
    print '<br>unfolded state 1 node number = <a href="isegment.cgi?iseg=n%s&dag=%s&">%s</a>'%(u1,dag,u1)
    if cuttype != "SEAM":
        print '<br>unfolded state 2 node number = <a href="isegment.cgi?iseg=n%s&dag=%s&">%s</a>'%(u2,dag,u2)
    if 'h' or 'p' or 'b' in ctype:
        print '<br>%s: entropy = %8.3f'%(cuttype,float(ntrp))
    print '<br>Relative traffic =%8.3f'%(float(traf))
    if cuttype == 'SEAM':
        try:
            readDag = open(dag,'r')
        except IOError:
            oneliner('<font color="#FF0000">BUG! DAG file open statement in tstatelineoutput</font>'+dag)
            return
        print '<br>Seam %5i limits='%(int(iseam))
        iseam = abs(int(iseam))
        for line in readDag:
            if line[0:5]=="SEAM ":
                line = line.split()
                i = line[1]
                nb = line[2]
                nrg = line[3]
                x1 = line[4]
                x2 = line[5]
                y1 = line[6]
                y2 = line[7]
                if int(i)==iseam:
                    print '%5i%5i%5i%5i'%(int(x1),int(x2),int(y1),int(y3))
        readDag.close()

    
def istateline(line,bline,dag):
    tics="....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|"
    nums="1--------10--------20--------30--------40--------50--------60--------70--------80--------90-------100"
    line = line.split()
    nseg = line[1]
    i = line[2]
    nsym = line[3]
    sas = line[4]
    ntrp = line[5]
    nvoid = line[6]
    nhb = line[7]
    ftype = line[17]
    conc = line[8]
    barrels = line[9:17]
    
    print '<font color="#11AA55"<pre>%s</pre></font>'%(" ".join(line))
    j = bline.find(" ")-1
    print '<font color="#112211"><pre>%s</pre></font>'%(nums)
    print '<font color="#112211"><pre>%s</pre></font>'%(tics)
    for step in range(0,j,100):
        nextstep = step + 100
        if nextstep > j:
            nextstep = j
        print '<font color="#AA1155"><pre>%s %s</pre></font>'%(bline[step:nextstep],nextstep)
    print '<table border=1>'
    print '<tr><td bgcolor="#BBAAFF">Segment node</td><td bgcolor="#FE8877">%9i</td></tr>'%(int(nseg))
    print '<tr><td bgcolor="#CC99FF">Folded state</td>'
    ftypes = {"F":'<td bgcolor="#CAFABA">FOLDED</td></tr>','I':'<td bgcolor="#BACAFA">INTERMEDIATE</td></tr>','U':'<td bgcolor="#FACADA">UNFOLDED</td></tr>'}
    print ftypes[ftype]
    nres = 0
    nch = 0
    lstch = '.'
    for step in range(0,j):
        if bline[step] != '.':
            nres += 1
            if lstch != bline[step]:
                nch += 1
        lstch = bline[step]
    print '<tr><td bgcolor="#BBAAFF">Size     </td><td bgcolor="#DACAFA"> %s </td></tr>'%(nres)
    print '<tr><td bgcolor="#AABBFF">Number of chains </td><td bgcolor="#CADAFA"> %s </td></tr>'%(nch)
    print '<tr><td bgcolor="#BBAAFF">Symmetry-generated copies    </td><td bgcolor="#DACAFA"> %s </td></tr>'%(nsym)
    print '<tr><td bgcolor="#AABBFF">Total buried surface area    </td><td bgcolor="#CADAFA"> %9.2f </td></tr>'%(float(sas))
    print '<tr><td bgcolor="#BBAAFF">Total unexpressed sidechain entropy    </td><td bgcolor="#DACAFA"> %9.2f </td></tr>'%(float(ntrp))
    print '<tr><td bgcolor="#AABBFF">Total remaining buried voids </td><td bgcolor="#CADAFA"> %s </td></tr>'%(nvoid)
    print '<tr><td bgcolor="#BBAAFF">Total remaining H-bonds </td><td bgcolor="#DACAFA"> %s </td></tr>'%(nhb)
    print '<tr><td bgcolor="#AABBFF">Concentration at end of simulation</td><td bgcolor="#CADAFA"> %8.2e </td></tr>'%(float(conc))
    print '<tr><td bgcolor="#BBAAFF">Barrels: </td><td bgcolor="#DACAFA">%s %s %s %s %s %s %s %s</td></tr>'%tuple(barrels)  
    print '</table>'
    try:
        readDag = open(dag,'r')
    except IOError:
        oneliner("Error opening dagfile line 137" + dag)
    print '<h5>Transition states (folding)</h5><pre>'
    print '       %7s%7s%7s%7s'%('tstate','f','u1','u2')
    for line in readDag:
        if line[0:6] == 'TSTATE':
            oldline = line
            line = line.split()
            nts = line[1]
            f = line[2]
            u1 = line[3]
            u2 = line[4]
            if u1 == nseg or u2 == nseg:
                outline = line[0]                
                outline += '<a href="isegment.cgi?iseg=tu%s&dag=%s&">%7i</a>'%(nts,dag,int(nts))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(f,dag,int(f))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(u1,dag,int(u1))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(u2,dag,int(u2))
                outline += " "+oldline[35:].strip()
                print outline
    print '</pre>'
    readDag.close()
    readDag = open(dag,'r')
    print '<h5>Transition state (unfolding)</h5><pre>'
    print '       %7s%7s%7s%7s'%('tstate','f','u1','u2')
    for line in readDag:
        if line[0:6] == 'TSTATE':
            oldline = line
            line = line.split()
            nts = line[1]
            f = line[2]
            u1 = line[3]
            u2 = line[4]
            if f == nseg:
                outline = line[0]                
                outline += '<a href="isegment.cgi?iseg=tu%s&dag=%s&">%7i</a>'%(nts,dag,int(nts))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(f,dag,int(f))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(u1,dag,int(u1))
                outline += '<a href="isegment.cgi?iseg=n%s&dag=%s&">%7i</a>'%(u2,dag,int(u2))
                outline += " "+oldline[35:].strip()
                print outline
    readDag.close()
    print '</pre>'
            

def oneliner(entry):
    print("<font color=\"#000055\"><pre>%s</pre></font>"%(entry))     
    
#HTML header
print "Content-Type: text/html;charset=utf-8"
print


print("<html><head>")
print("</head>")
query = os.environ['QUERY_STRING']
iseg = parseit(query,"iseg")
dag = parseit(query,"dag")
print '<br>'
readDagFile(iseg,dag)
print("</body></html>")

