
from math import *
import random



def roll( n_rolls=0 ):
    
    return [ random.randint(1,6) for k in range(n_rolls) ]



def dumb_defend( n_defending, attack_rolls ):
    
    if n_defending == 1:
        return 1
    
    if n_defending >= 2:
        return 2
    
    return 0

def dumb_attack( n_attacking ):
    
    if n_attacking == 2:
        return 1
    
    if n_attacking == 3:
        return 2
    
    if n_attacking >= 4:
        return 3
    
    return 0

def once_outcome( n_attacking, n_defending, attack_strategy, defend_strategy ):
    
    attack_rolls = roll(dumb_attack( n_attacking ))
    defend_rolls = roll(dumb_defend(n_defending, attack_rolls))
    
    attack_rolls.sort()
    defend_rolls.sort()
    attack_rolls.reverse()
    defend_rolls.reverse()
    
    attack_loss = 0
    defend_loss = 0
    
    if defend_rolls == [] or attack_rolls == []:
        return (0,0)
    
    if defend_rolls[0] >= attack_rolls[0]:
        attack_loss += 1
    else:
        defend_loss += 1
    
    if len(defend_rolls) > 1 and len(attack_rolls) > 1:
        if defend_rolls[1] >= attack_rolls[1]:
            attack_loss += 1
        else:
            defend_loss += 1
    
    return ( defend_loss, attack_loss )



def outcome( n_attacking, n_defending, attack_strategy=dumb_attack, defend_strategy=dumb_defend ):
    
    attacking = n_attacking
    defending = n_defending
    
    while attacking>1 and defending>0:
        
        ( defend_loss, attack_loss ) = once_outcome( attacking, defending, attack_strategy, defend_strategy )
        
        #print (defending, attacking ), ( defend_loss, attack_loss )
        
        attacking -= attack_loss
        defending -= defend_loss
    
    return ( defending, attacking )



def atackerwin_probability( n_attacking, n_defending, attack_strategy=dumb_attack, defend_strategy=dumb_defend, n_trials=1000 ):
    
    attack_wins = 0
    
    for n in range(n_trials):
        
        ( defending, attacking ) = outcome( n_attacking, n_defending, attack_strategy, defend_strategy )
        
        if defending == 0:
            attack_wins += 1
    
    return float(attack_wins)/n_trials



if __name__ == "__main__":
    
    
    for k in range(5):
        print atackerwin_probability( 10,20,  n_trials = 10000)
    
    










