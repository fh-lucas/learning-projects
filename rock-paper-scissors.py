import random

rock = '''
    _______
---'   ____)
      (_____)
      (_____)
      (____)
---.__(___)
'''

paper = '''
    _______
---'   ____)____
          ______)
          _______)
         _______)
---.__________)
'''

scissors = '''
    _______
---'   ____)____
          ______)
       __________)
      (____)
---.__(___)
'''

game = [rock, paper, scissors]
player = int(input("Choose 0 for Rock, 1 for Paper or 3 for Scissors\n"))
print("\nPlayer:" + game[player])
machine = random.randint(0, len(game)-1)
print("Machine:" + game[machine])

if player == 0:
    if machine == 0:
        print("Tie")
    elif machine == 1:
        print("You lost")
    else:
        print("You win")
elif player == 1:
    if machine == 0:
        print("You win")
    elif machine == 1:
        print("Tie")
    else:
        print("You lost")
else:
    if machine == 0:
        print("You lost")
    elif machine == 1:
        print("You win")
    else:
        print("Tie")