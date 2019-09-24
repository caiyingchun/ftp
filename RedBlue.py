from time import time
from random import seed, randint
from tkinter import Tk, Label, Entry, Button
seed(time())

def yaner():
    for num_r_ball in nums_r_ball:
        num_r_ball.delete(0, 'end')
        num_r_ball.insert(0, randint(1, 33))
    num_b_ball.delete(0, 'end')
    num_b_ball.insert(0, randint(1, 16))
    
root = Tk()
root.title('小c加油！！！')
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
win_width, win_height = 415, 210
x = (screen_width - win_width) / 2
y = (screen_height - win_height) / 2
root.geometry('%dx%d+%d+%d' %(win_width, win_height, x, y))
root.resizable(0,0)

lb_title = Label(root, text='双色球机选器', width=15, bg='lightgrey', fg='purple', font=('Arial', 24))
lb_title.grid(row=0, column=0, columnspan=7, padx=10, pady=20, sticky='WENS')

lb_red = Label(root, text='红球', fg='red', font=('Arial', 15))
lb_red.grid(row=1, column=0, sticky='W')

lb_blue = Label(root, text='蓝球', fg='blue', font=('Arial', 15))
lb_blue.grid(row=1, column=6, sticky='W')

nums_r_ball = []
for idx in range(6):
    num_r_ball = Entry(root, width=5, fg='red', font=('Arial', 12), justify='center')
    num_r_ball.grid(row=2, column=idx, padx=5, pady=3)
    nums_r_ball.append(num_r_ball)

num_b_ball = Entry(root, width=5, fg='blue', font=('Arial', 12), justify='center')
num_b_ball.grid(row=2, column=6, padx=5, pady=3)

bt = Button(root, text='摇号', command=yaner, width=10, activeforeground='Magenta', fg='Magenta', font=('Arial', 18), relief='groove')
bt.grid(row=3, column=1, columnspan=5, pady=10, ipadx=15)

root.mainloop()
