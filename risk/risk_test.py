
from risk import *




defend_map = { 
(1,1):2,
(2,1):2,
(3,1):2,
(4,1):2,
(5,1):2,
(6,1):2,

(2,2):2,
(3,2):2,
(4,2):2,
(5,2):2,
(6,2):2,

(3,3):2,
(4,3):2,
(5,3):2,
(6,3):2,

(4,4):1,
(5,4):1,
(6,4):1,

(5,5):1,
(6,5):1,

(6,6):1 } 



def defend( n_defending, attack_rolls ):
    
    print attack_rolls
    if n_defending == 1:
        return 1
    
    if n_defending >= 2:
            
        atack_rolls.sort()
        attack_rolls.reverse()
        
        #return defend_map( (attack_rolls[0],attack_rolls[1]) )
    
    
    print 'problem'
    return 0




if __name__ == "__main__":
    
    
    for k in range(5):
        print atackerwin_probability( n_attacking=4,n_defending=4,  n_trials = 1000, defend_strategy =defend )
    
    








