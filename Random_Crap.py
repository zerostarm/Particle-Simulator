'''
Created on Jul 5, 2022

@author: Stephen
'''

'''
            r_small = True
            with np.errstate(divide='raise', invalid='raise'):
                try:
                    r =  r*r_min*2/r_norm
                except:
                    r = np.asarray([1,0,0])*r_min*2
            temp_position = particle2.getPosition() + r
            r = temp_position - particle2.getPosition(u)
            r_norm = r_min**2
            r_hat = r/r_norm'''
#if particle2.Index < particle1.Index:
        #    u = -2
        #r = particle1.getPosition() - particle2.getPosition(u)
  
