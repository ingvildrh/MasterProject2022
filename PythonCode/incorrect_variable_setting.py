if (i2 <= 0):
                i2= []
            if (i2): 
                iz = i2z
                actset[iz] = 1 #hvorfor setter denne her to variabler, både actset og actset0??
                qc = -1
            else:
                iz = []
        
        if (iz):
            qu = IGIs[:, iz]
            vA = np.transpose(np.matrix(Qmat0i@(qu)))
            qdiv = qc+vA[iz] 

        ...
    ...
