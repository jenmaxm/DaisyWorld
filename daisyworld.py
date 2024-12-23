import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import *
import tkinter.font as tkFont
from random import randint
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

class main:
  def __init__(self):
    self.count = 0
    self.S = 9.17 * (10**5)
    self.s = 5.75 * (10**-5)
    self.g = 0.3
    self.Aw = 0.75
    self.Ab = 0.25
    self.q = 2.06 * (10**9)
    self.max_L = 1.60
    self.L = np.arange(0.6,self.max_L,0.001)
    self.Ti_list = np.arange(280,310,0.03)

  def Te_list(self):
    self.Te_list = []
    for i in range(len(self.L)):
      Te =  pow((self.S*self.L[i])/(2*self.s), 1/4)
      self.Te_list.append(Te)

  def growth_rate(self):
    self.b_list = []
    for i in range(len(self.Ti_list)):
      b = 1 - (0.003265 * ((295.5 - self.Ti_list[i]) ** 2 ))
      self.b_list.append(b)

  def tb_method(self):
    Tb_list = []
    for i in range(len(self.L)):
      ab = (((17.5 ** 2) * (1-self.g)) - (((self.Ti_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Ti_list[i]-295.5) ** 2 ))
      if ab >= 0 and ab <= 0.7:
        Tb = pow((0.25 * self.q) + ((0.5 * self.S * self.L[i]) / (self.s)) + (ab * ((0.25*self.S*self.L[i]) / self.s) - (0.25 * self.q)), 1/4)
        if Tb > 270 and Tb < 320:
          Tb_list.append(Tb)
    self.Tb_list = Tb_list

  def luminosity_black_method(self):
    self.tb_method()
    self.Lb_list = []
    self.ab_list = []
    for i in range(len(self.Tb_list)):
      ab = (((17.5 ** 2) * (1-self.g)) - (((self.Tb_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tb_list[i]-295.5) ** 2 ))
      if ab >= 0 and ab <= 0.7:
        Lb = ((self.Tb_list[i] ** 4) + (0.25* self.q * ab) - (0.25*self.q)) / (((0.5*self.S)/self.s)) + ((ab) * ((0.25*self.S)/self.s))
        self.Lb_list.append(Lb)
        self.ab_list.append(ab)

  def mean_temp_black(self):
    self.tb_method()
    self.Te_b_list = []
    for i in range(len(self.Tb_list)):
      self.ab = (((17.5 ** 2) * (1-self.g)) - (((self.Tb_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tb_list[i]-295.5) ** 2 ))
      if self.ab >= 0 and self.ab <= 0.7:
        Te_b = pow(((self.S*self.L[i])/(self.s) * (0.5 + (0.25*self.ab))), 1/4 )
        self.Te_b_list.append(Te_b)
      else:
        self.ab = 0
        Te_b = pow(((self.S*self.L[i])/(self.s) * (0.5 + (0.25*self.ab))), 1/4 )
        self.Te_b_list.append(Te_b)

  def luminosity_white_method(self):
    self.tw_method()
    self.Lw_list = []
    self.aw_list = []
    for i in range(len(self.Tw_list)):
      aw = (((17.5 ** 2) * (1-self.g)) - (((self.Tw_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tw_list[i]-295.5) ** 2 ))
      if aw >= 0 and aw <= 0.7:
        Lw = ((self.Tw_list[i] ** 4) - (0.25* self.q * aw) + (0.25*self.q)) / (((0.5*self.S)/self.s)) - ((aw) * ((0.25*self.S)/self.s))
        self.Lw_list.append(Lw)
        self.aw_list.append(aw)

  def mean_temp_white(self):
    self.tw_method()
    self.Te_w_list = []
    for i in range(len(self.Tw_list)):
      self.aw = (((17.5 ** 2) * (1-self.g)) - (((self.Tw_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tw_list[i]-295.5) ** 2 ))
      if self.aw >= 0 and self.aw <= 0.7:
        Te_w = pow(((self.S*self.L[i])/(self.s) * (0.5 - (0.25*self.aw))), 1/4 )
        self.Te_w_list.append(Te_w)
      else:
        self.aw = 0
        Te_w = pow(((self.S*self.L[i])/(self.s) * (0.5 - (0.25*self.aw))), 1/4 )
        self.Te_w_list.append(Te_w)
    
  def two_species_temperature(self):
    self.Te_b2_list = []
    self.ab_2_list = []
    for i in range(len(self.L)):
      ab_2 = (0.829 / (self.L[i] - 0.129)) - 0.664
      if ab_2 >=0 and ab_2 <= 0.673: 
        Te_b2 = pow((self.S*self.L[i])/(self.s) * (0.332 + ab_2), 1/4 )
        self.Te_b2_list.append(Te_b2)
        self.ab_2_list.append(ab_2)
      else:
        ab_2 = 0
        Te_b2 = pow((self.S*self.L[i])/(self.s) * (0.332 + ab_2), 1/4 )
        self.Te_b2_list.append(Te_b2)
        self.ab_2_list.append(ab_2)
    self.ab_list = self.ab_2_list
    self.sorting_ab()

  def two_species_area(self):
    self.Te_b2_list = []
    self.ab_2_list = []
    for i in range(len(self.L)):
      ab_2 = (0.829 / (self.L[i] - 0.129)) - 0.664
      if ab_2 >=0 and ab_2 <= 0.673: 
        Te_b2 = pow(((self.S*self.L[i])/(self.s) * (0.5 + (0.25*ab_2))), 1/4 )
        self.Te_b2_list.append(Te_b2)
        self.ab_2_list.append(ab_2)
      else:
        ab_2 = 0
        Te_b2 = pow(((self.S*self.L[i])/(self.s) * (0.5 + (0.25*ab_2))), 1/4 )
        self.Te_b2_list.append(Te_b2)
        self.ab_2_list.append(ab_2)
    self.ab_list = self.ab_2_list
    self.sorting_ab()

    self.Te_w2_list = []
    self.aw_2_list = []
    for i in range(len(self.L)):
      aw_2 = 0.7 - self.ab_2_list[i] 
      if aw_2 >=0 and aw_2 <= 0.673:
        Te_w2 = pow(((self.S*self.L[i])/(self.s) * (0.5 - (0.25*aw_2))), 1/4 )
        self.Te_w2_list.append(Te_w2)
        self.aw_2_list.append(aw_2)
      else:
        aw_2 = 0
        Te_w2 = pow(((self.S*self.L[i])/(self.s) * (0.5 - (0.25*aw_2))), 1/4 )
        self.Te_w2_list.append(Te_w2)
        self.aw_2_list.append(aw_2)
    self.aw_list = self.aw_2_list
    self.sorting_aw()

  def tw_method(self):
    Tw_list = []
    for i in range(len(self.L)):
      aw = (((17.5 ** 2) * (1-self.g)) - (((self.Ti_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Ti_list[i]-295.5) ** 2 ))
      if aw >= 0 and aw <= 0.7:
        Tw = pow((-0.25 * self.q) + ((0.5 * self.S * self.L[i]) / (self.s)) + (aw * ((-1*0.25*self.S*self.L[i]) / self.s) + (0.25 * self.q)), 1/4)
        if Tw > 270 and Tw < 320:
          Tw_list.append(Tw)
    self.Tw_list = Tw_list
    print(self.Tw_list)
      
  def white_area(self):
    self.tw_method()
    Lw_list = []
    aw_list = []
    for i in range(len(self.Tw_list)):
      aw = (((17.5 ** 2) * (1-self.g)) - (((self.Tw_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tw_list[i]-295.5) ** 2 ))
      if aw >= 0 and aw <= 0.7:
        Lw = ((self.Tw_list[i] ** 4) - (0.25* self.q * aw) + (0.25*self.q)) / (((0.5*self.S)/self.s)) - ((aw) * ((0.25*self.S)/self.s))
        Lw_list.append(Lw)
        aw_list.append(aw)
    self.aw_list = aw_list
    self.sorting_aw()

  def black_area(self):
    self.tb_method()
    Lb_list = []
    ab_list = []
    for i in range(len(self.Tb_list)):
      ab = (((17.5 ** 2) * (1-self.g)) - (((self.Tb_list[i]-295.5) ** 2 ))) / ((17.5 ** 2) - ((self.Tb_list[i]-295.5) ** 2 ))
      if ab >= 0 and ab <= 0.7:
        Lb = ((self.Tb_list[i] ** 4) + (0.25* self.q * ab) - (0.25*self.q)) / (((0.5*self.S)/self.s)) + ((ab) * ((0.25*self.S)/self.s))
        Lb_list.append(Lb)
        ab_list.append(ab)
    self.ab_list = ab_list
    
    self.sorting_ab()

  def sorting_ab(self):
    x = self.ab_list
    a = 200
    gui_ab = []
    b = 0
    for i in range(a):
      try:
        gui_ab.append(x[b])
        b = b + 4
      except:
          pass
    for i in range(len(gui_ab)):
        gui_ab[i] = gui_ab[i] * 100
        gui_ab[i] = int(gui_ab[i])
    self.gui_ab = gui_ab

  def sorting_aw(self):
    x = self.aw_list
    a = 200
    gui_aw = []
    b = 0
    for i in range(a):
      try:
        gui_aw.append(x[b])
        b = b + 4
      except:
          pass
    for i in range(len(gui_aw)):
        gui_aw[i] = gui_aw[i] * 100
        gui_aw[i] = int(gui_aw[i])
    self.gui_aw = gui_aw

  def mul_plots(self,x,y,y2,xl,yl):
    root = Tk()
    root.title('Graph')
    root.geometry('700x700')

    fig = Figure(figsize = (5, 5),dpi = 100)
    ax1 = Figure(figsize = (5,5), dpi= 100)
    
    fig, ax1 = plt.subplots()
    ax1.plot(x,y,color='blue')
    ax2 = ax1.twinx()
    ax2.plot(x,y2,color='green')
    ax1.set_xlabel(xl)
    ax1.set_ylabel(yl)

    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()

  def plot(self,x,y,xl,yl):
    root = Tk()
    root.title('Graph')
    root.geometry('500x500')

    fig = Figure(figsize = (5, 5),dpi = 100)
    
    plot1 = fig.add_subplot()
    plot1.plot(x,y)
    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()
    plot1.set_xlabel(xl)
    plot1.set_ylabel(yl)

    root.mainloop()

  def clearFrame(self):
    for widget in self.container.winfo_children():
       widget.destroy()

    self.container.configure(bg="saddle brown")
    self.container.grid_forget()

class GUI(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Main Menu")

    times36 = tkFont.Font(family="Times New Roman",size=36,weight="bold")
    times24 = tkFont.Font(family="Times New Roman",size=18,weight="bold")

    self.container = Frame(master,width=500,height=500,bg='grey')
    self.container.grid()

    self.label = Label(self.container, text='Daisyworld Model',font=times36,bg='grey')
    self.label.grid(row=0,column=0,pady=5,padx=10)

    self.bbutton = Button(self.container, text="Black Species Only",command=self.black,font=times24,relief=RAISED)
    self.bbutton.grid(row = 1,column = 0,pady=5)

    self.wbutton = Button(self.container, text="White Species Only",command=self.white,font=times24,relief=RAISED)
    self.wbutton.grid(row = 2,column = 0,pady=5)

    self.bwbutton = Button(self.container, text="Black and White Species",command=self.black_white,font=times24,relief=RAISED)
    self.bwbutton.grid(row = 3,column = 0,pady=5)

    self.fpbutton = Button(self.container, text="Other Graphs",command=self.graphs,font=times24,relief=RAISED)
    self.fpbutton.grid(row = 4,column = 0,pady=5)

    self.cvbutton = Button(self.container, text="Change Starting Values",command=self.change_values,font=times24,relief=RAISED)
    self.cvbutton.grid(row = 5,column = 0,pady=5)

  def black(self):
      self.clearFrame()
      black_GUI(root)

  def white(self):
      self.clearFrame()
      white_GUI(root)

  def black_white(self):
      self.clearFrame()
      two_speciesGUI(root)

  def graphs(self):
      self.clearFrame()
      other_plots(root)

  def change_values(self):
      self.clearFrame()
      change_values_GUI(root)
 
class black_GUI(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Area of Black Daisies(Single species)")

    self.times18 = tkFont.Font(family="Times New Roman",size=18,weight="bold")

    self.container = Frame(master,bg='grey')
    self.container.grid()

    self.canvas = Canvas(self.container,width=400,height=400,bg='saddle brown')
    self.canvas.grid()
       
    self.show_button = Button(self.container, text="Show timestep",command=self.b_area,font=self.times18)
    self.show_button.grid(row = 1,column = 0,pady=5)

    self.skip_button = Button(self.container, text="Next timestep", font=self.times18,command=self.clicked)
    self.skip_button.grid(row = 2 , column = 0,pady=5)

    self.skip_button_10 = Button(self.container, text="Skip 10 timesteps", command=self.clicked_10,font=self.times18)
    self.skip_button_10.grid(row = 3 , column =0,pady=5)

    self.button = Button(self.container, text="Explain Please", command=self.explain,font=self.times18)
    self.button.grid(row =4,column =0,pady=5)

    self.plot_button = Button(self.container,text = "Plot",command = self.show_plot,font=self.times18)
    self.plot_button.grid(row=5,column=0,pady=5)

    self.return_button = Button(self.container,text='Exit',command = self.exit,font=self.times18)
    self.return_button.grid(row=6,column=0,pady=5)

    self.label = Label(self.container,width = 40,height=3,text=None,bg='grey',font=self.times18)
    self.label.grid(row=7,column= 0)
  

  def plot(self,x,y,xl,yl):

    root = Tk()
    root.title('Graph')
    root.geometry('500x500')

    fig = Figure(figsize = (5, 5),dpi = 100)
    
    plot1 = fig.add_subplot(111)
    plot1.plot(x,y)
    plot1.set_xlabel(xl)
    plot1.set_ylabel(yl)
    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()

    root.mainloop()

  def show_plot(self):
    self.tb_method()
    self.black_area()
    self.plot(self.L[:len(self.ab_list)],self.ab_list,'Luminosity','Area of Black Daisies')

  def black_daisy(self):
    """Creates a single black daisy (circle) on the canvas."""
    x = randint(0, 350)
    y = randint(0, 350)
    r = 5
    self.canvas.create_oval(x - r, y - r, x + r, y + r, fill='black')

  def b_area(self):
    """Updates the canvas to display black daisies covering the specified percentage of the area."""
    self.black_area()  # Assuming this is defined elsewhere to handle setup tasks.
    self.canvas.delete("all")

    if self.count > (len(self.gui_ab) - 1):
        self.count = 0

    else:
        grid_area = 350 * 350  # Total grid area.
        percentage = self.gui_ab[self.count] / 100  # Convert percentage to a decimal.
        total_covered_area = grid_area * percentage  # Calculate the area to be covered.

        daisy_area = np.pi * (5 ** 2)  # Area of one daisy (circle with radius 15).
        number_of_daisies = int(total_covered_area / daisy_area)  # Calculate number of daisies needed.

        for i in range(number_of_daisies):
            self.black_daisy()

        percent = f"{self.gui_ab[self.count]}% Area Covered by Black Daisies"
        self.label.configure(text=percent)

    self.show_button["state"] = "disabled"

  
  def clicked(self):
    self.count 
    self.count = self.count + 1
    self.show_button["state"] = NORMAL

  def clicked_10(self):
    self.count 
    self.count = self.count + 10
    self.show_button["state"] = NORMAL

  def explain(self):
    self.configfile = Text(self.container,width=80,bg='black')
    self.configfile.grid(row=0,column=2)
    file = open('black_exp.txt')
    self.configfile.insert(INSERT, file.read())
    self.configfile.config(state= DISABLED)

    self.clear_button = Button(self.container,command=self.clear_label,text='Clear Explanation',font=self.times18)
    self.clear_button.grid(row=1,column=2)

  def clear_label(self):
    self.configfile.after(100, self.configfile.destroy())
    self.clear_button.after(100, self.clear_button.destroy)

  def exit(self):
      self.clearFrame()
      GUI(root)
    

class white_GUI(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Area of White Daisies(Single species)")

    self.times18 = tkFont.Font(family="Times New Roman",size=18,weight="bold")

    self.container = Frame(master, bg='grey')
    self.container.grid()

    self.canvas = Canvas(self.container,width=400,height=400,bg='saddle brown')
    self.canvas.grid()
       
    self.show_button = Button(self.container, text="Show timestep", command=self.w_area,font=self.times18)
    self.show_button.grid(row = 1,column = 0,pady=5)

    self.skip_button = Button(self.container, text="Next timestep", command=self.clicked, font=self.times18)
    self.skip_button.grid(row = 2 , column = 0,pady=5)

    self.skip_button_10 = Button(self.container, text="Skip 10 timesteps", command=self.clicked_10, font=self.times18)
    self.skip_button_10.grid(row = 3 , column =0,pady=5)

    self.button = Button(self.container, text="Explain Please", command=self.explain, font=self.times18)
    self.button.grid(row =4,column =0,pady=5)

    plot_button = Button(self.container,text = "Plot",command = self.show_plot, font=self.times18)
    plot_button.grid(row=5,column=0,pady=5)
    
    self.return_button = Button(self.container,text='Exit',command = self.exit,font=self.times18)
    self.return_button.grid(row=6,column=0,pady=5)

    self.label = Label(self.container,width = 40,height=3,text=None,bg='grey',font=self.times18)
    self.label.grid(row=7,column= 0,pady=5)
  

  def plot(self,x,y,xl,yl):

    root = Tk()
    root.title('Graph')
    root.geometry('500x500')

    fig = Figure(figsize = (5, 5),dpi = 100)
    
    plot1 = fig.add_subplot()
    plot1.plot(x,y)
    plot1.set_xlabel(xl)
    plot1.set_ylabel(yl)
    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()

    root.mainloop()

  def show_plot(self):
    self.tw_method()
    self.white_area()
    self.plot(self.L[:len(self.aw_list)],self.aw_list,'Luminosity','Area of White Daisies')

  def white_daisy(self):
    """Creates a single black daisy (circle) on the canvas."""
    x = randint(0, 350)
    y = randint(0, 350)
    r = 5
    self.canvas.create_oval(x - r, y - r, x + r, y + r, fill='white')

  def w_area(self):
    """Updates the canvas to display white daisies covering the specified percentage of the area."""
    self.white_area() ]
    self.canvas.delete("all")

    if self.count > (len(self.gui_aw) - 1):
        self.count = 0

    else:
        grid_area = 350 * 350 
        percentage = self.gui_aw[self.count] / 100 
        total_covered_area = grid_area * percentage 

        daisy_area = np.pi * (5 ** 2)
        number_of_daisies = int(total_covered_area / daisy_area)  

        for i in range(number_of_daisies):
            self.white_daisy()

        percent = f"{self.gui_aw[self.count]}% Area Covered by White Daisies"
        self.label.configure(text=percent)

    self.show_button["state"] = "disabled"
  
  def clicked(self):
    self.count 
    self.count = self.count + 1
    self.show_button["state"] = NORMAL

  def clicked_10(self):
    self.count 
    self.count = self.count + 10
    self.show_button["state"] = NORMAL

  def explain(self):
    self.configfile = Text(self.container,width=100,bg='black')
    self.configfile.grid(row=0,column=2)
    file = open('white_exp.txt')
    self.configfile.insert(INSERT, file.read())
    self.configfile.config(state= DISABLED)

    self.clear_button = Button(self.container,command=self.clear_label,text='Clear Explanation',font=self.times18)
    self.clear_button.grid(row=1,column=2)

  def clear_label(self):
    self.configfile.after(100, self.configfile.destroy())
    self.clear_button.after(100, self.clear_button.destroy)

  def exit(self):
      self.clearFrame()
      GUI(root)


class two_speciesGUI(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Area of Both Daisies")

    self.times18 = tkFont.Font(family="Times New Roman",size=18,weight="bold")

    self.container = Frame(master, bg='grey')
    self.container.grid()

    self.canvas = Canvas(self.container,width=400,height=400,bg='saddle brown')
    self.canvas.grid()
       
    self.show_button = Button(self.container, text="Show timestep", command=self.b_w_area,font=self.times18)
    self.show_button.grid(row = 1,column=0,pady=5)

    self.skip_button = Button(self.container, text="Next timestep", command=self.clicked,font=self.times18)
    self.skip_button.grid(row = 2 ,column= 0,pady=5)

    self.skip_button = Button(self.container, text="Skip 10 timesteps", command=self.clicked_10,font=self.times18)
    self.skip_button.grid(row =3,column =0,pady=5)
    
    self.return_button = Button(self.container,text='Exit',command = self.exit,font=self.times18)
    self.return_button.grid(row=6,column=0,pady=5)

    self.label = Label(self.container,width = 40,height=3,text=None,bg='grey',font=self.times18)
    self.label.grid(row=7,column= 0,pady=5)

    self.label2 = Label(self.container,width = 40,height=3,text=None,bg='grey',font=self.times18)
    self.label2.grid(row=8,column= 0,pady=5)

    self.button = Button(self.container, text="Plot", command=self.show_plot,font=self.times18)
    self.button.grid(row =5,column =0,pady=5)

    self.button = Button(self.container, text="Explain Please", command=self.explain,font=self.times18)
    self.button.grid(row =4,column =0,pady=5)

  def mul_plots(self,x,y,y2,yl,xl):
    root = Tk()
    root.title('Graph')
    root.geometry('700x700')

    fig = Figure(figsize = (5, 5),dpi = 100)
    
    fig, ax1 = plt.subplots()
    ax1.plot(x,y,color='blue')
    ax2 = ax1.twinx()
    ax2.plot(x,y2,color='green')
    ax1.set_xlabel(xl)
    ax1.set_ylabel(yl)

    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()

  def show_plot(self):
    self.two_species_area()
    self.mul_plots(self.L[:len(self.ab_list)],self.ab_list,self.aw_list,'Total Area','Luminosity')

  def clicked(self):
    self.count 
    self.count = self.count + 1
    self.show_button["state"] = NORMAL

  def clicked_10(self):
    self.count 
    self.count = self.count + 10
    self.show_button["state"] = NORMAL

  def white_daisy(self):
    x = randint(0,350)
    y = randint(0,350)
    r = 5
    self.dot = self.canvas.create_oval(x-r,y-r,x+r,y+r,fill='white')
    return self.dot

  def black_daisy(self):
    x = randint(0,350)
    y = randint(0,350)
    r = 5
    self.dot = self.canvas.create_oval(x-r,y-r,x+r,y+r,fill='black')
    return self.dot

  def b_w_area(self):
      """Updates the canvas to display black and white daisies covering specified percentages of the area."""
      self.two_species_area() 
      self.canvas.delete("all")

      if self.count > (len(self.gui_ab) - 1) or self.count > (len(self.gui_aw) - 1):
          self.count = 0

      else:
          grid_area = 350 * 350  # Total grid area.

          # Black daisies
          percentage_b = self.gui_ab[self.count] / 100
          total_covered_area_b = grid_area * percentage_b
          daisy_area = np.pi * (5 ** 2)
          number_of_black_daisies = int(total_covered_area_b / daisy_area)

          for i in range(number_of_black_daisies):
              self.black_daisy()

          # White daisies
          percentage_w = self.gui_aw[self.count] / 100
          total_covered_area_w = grid_area * percentage_w
          number_of_white_daisies = int(total_covered_area_w / daisy_area)

          for i in range(number_of_white_daisies):
              self.white_daisy()

          percent_b = f"{self.gui_ab[self.count]}% : of Area Covered in Black Daisies"
          percent_w = f"{self.gui_aw[self.count]}% : of Area Covered in White Daisies"
          self.label.configure(text=percent_b)
          self.label2.configure(text=percent_w)

      self.show_button["state"] = "disabled"

  def explain(self):
    self.configfile = Text(self.container,width=100)
    self.configfile.grid(row=0,column=2)
    file = open('both_exp.txt')
    self.configfile.insert(INSERT, file.read())
    self.configfile.config(state= DISABLED)

    self.clear_button = Button(self.container,command=self.clear_label,text='Clear Explanation',font=self.times18)
    self.clear_button.grid(row=1,column=2)

  def clear_label(self):
    self.configfile.after(100, self.configfile.destroy())
    self.clear_button.after(100, self.clear_button.destroy)

  def exit(self):
    self.clearFrame()
    GUI(root)
    

class other_plots(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Graphs")

    self.times18 = tkFont.Font(family="Times New Roman",size=18,weight="bold")

    self.container = Frame(master,bg='grey')
    self.container.grid()

    self.button = Button(self.container, text="Temperature of lifeless planet as luminosity increases", command=self.show_plot1,font=self.times18)
    self.button.grid(row =1,column =0,pady=5)

    self.button = Button(self.container, text="Area occupied by single daisies species against local temperature", command=self.show_plot3,font=self.times18)
    self.button.grid(row =2,column =0,pady=5)

    self.button = Button(self.container, text="Growth Rate against local temperature", command=self.show_plot2,font=self.times18)
    self.button.grid(row =3,column =0,pady=5)

    self.button = Button(self.container, text="Local temperature for black daisies against solar luminosity", command=self.show_plot4,font=self.times18)
    self.button.grid(row =4,column =0,pady=5)

    self.button = Button(self.container, text="Solar Luminonisity against local temperature for black daisies only species", command=self.show_plot5,font=self.times18)
    self.button.grid(row =5,column =0,pady=5)

    self.button = Button(self.container, text="Area of Black Daisises in one species model", command=self.show_plot6,font=self.times18)
    self.button.grid(row =6,column =0,pady=5)

    self.button = Button(self.container, text="Mean planetary temperature against solar luminosity for black daisies", command=self.show_plot7,font=self.times18)
    self.button.grid(row =7,column =0,pady=5)

    self.button = Button(self.container, text="Local temperature for white daisies against solar luminosity", command=self.show_plot8,font=self.times18)
    self.button.grid(row =8,column =0,pady=5)

    self.button = Button(self.container, text="Mean planetary temperature against solar luminosity for white daisies", command=self.show_plot9,font=self.times18)
    self.button.grid(row =9,column =0,pady=5)

    self.button = Button(self.container, text="Area of White Daisises in one species model", command=self.show_plot10,font=self.times18)
    self.button.grid(row =10,column =0,pady=5)

    self.button = Button(self.container, text="Black and White regulation", command=self.show_plot11,font=self.times18)
    self.button.grid(row =11,column =0,pady=5)

    self.return_button = Button(self.container,text='Exit',command = self.exit,font=self.times18)
    self.return_button.grid(row=12,column=0,pady=5)

  def mul_plots(self,x,y,y2,filename,xl,yl):
    root = Tk()
    root.title('Graph')
    root.geometry('700x700')

    fig = Figure(figsize = (5, 5),dpi = 100)
    ax1 = Figure(figsize = (5,5), dpi= 100)
    
    fig, ax1 = plt.subplots()
    ax1.plot(x,y,color='blue')
    ax2 = ax1.twinx()
    ax2.plot(x,y2,color='green')
    ax1.set_xlabel(xl)
    ax1.set_ylabel(yl)

    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()

    configfile = Text(root,width=100)
    configfile.pack()
    file = open(filename)
    configfile.insert(INSERT, file.read())
    configfile.config(state= DISABLED)

  def plot(self,x,y,filename,xl,yl):
    root = Tk()
    root.title('Graph')
    root.geometry('500x500')

    fig = Figure(figsize = (5, 5),dpi = 100)
    
    plot1 = fig.add_subplot()
    plot1.plot(x,y)
    canvas = FigureCanvasTkAgg(fig,master = root)  
    canvas.draw()
    canvas.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas,root)
    toolbar.update()
    canvas.get_tk_widget().pack()
    plot1.set_xlabel(xl)
    plot1.set_ylabel(yl)

    configfile = Text(root,width=100)
    configfile.pack()
    file = open(filename)
    configfile.insert(INSERT, file.read())
    configfile.config(state= DISABLED)

    root.mainloop()


  def show_plot1(self):
    self.Te_list()
    self.plot(self.L,self.Te_list,'plot1.txt','Luminosity','Mean Temperature of Planet')

  def show_plot2(self):
    self.growth_rate()
    self.plot(self.Ti_list[:len(self.b_list)],self.b_list,'plot2.txt','Temperature of Daisies','Growth Rate')

  def show_plot3(self): 
    self.growth_rate()
    a_list = []
    for i in range(len(self.Ti_list)):
      a = 1 - (self.g / self.b_list[i])
      if a >= 0 and a <= 0.7:
        a_list.append(a)
    self.plot(self.Ti_list[:len(a_list)], a_list,'plot3.txt','Temperature of Daisies','Area')

  def show_plot4(self):
    self.tb_method()
    self.plot(self.L[:len(self.Tb_list)],self.Tb_list,'plot4.txt','Luminosity','Temperature of Black Daisies')
    

  def show_plot5(self):
    self.luminosity_black_method()
    self.plot(self.Tb_list[:len(self.Lb_list)],self.Lb_list,'plot5.txt','Temperature of Black Daisies','Luminosity')


  def show_plot6(self):
    self.luminosity_black_method()
    self.plot(self.L[:len(self.ab_list)],self.ab_list,'plot6.txt','Luminosity','Area of black daisies')
  

  def show_plot7(self):
    self.mean_temp_black()
    self.plot(self.L[:len(self.Te_b_list)],self.Te_b_list,'plot7.txt','Luminosity','Mean Temperature of planet with black daisies')


  def show_plot8(self):
    self.luminosity_white_method()
    self.plot(self.Tw_list[:len(self.Lw_list)],self.Lw_list,'plot8.txt','Temperature of white daisies','Luminosity')


  def show_plot10(self):
    self.luminosity_white_method()
    self.plot(self.L[:len(self.aw_list)],self.aw_list,'plot9.txt','Luminosity','Area of white daisies')


  def show_plot9(self):
    self.mean_temp_white()
    self.plot(self.L[:len(self.Te_w_list)],self.Te_w_list,'plot10.txt','Luminosity','Mean Temperature of planet with white daisies')


  def show_plot11(self):
    self.two_species_temperature()
    self.plot(self.L[:len(self.Te_b2_list)],self.Te_b2_list,'plot11.txt','Luminosity','Temperature of planet')

  def exit(self):
      self.clearFrame()
      GUI(root)

class change_values_GUI(main):
  def __init__(self, master):
    super().__init__()
    self.master = master
    master.title("Change Values")

    self.times18 = tkFont.Font(family="Times New Roman",size=15,weight="bold")

    self.container = Frame(master)
    self.container.grid()
    
    self.death_rate = Scale(self.container, from_=0, to=1,digits = 3, resolution = 0.01, orient=HORIZONTAL,label = 'Death Rate',font=self.times18,width=40)
    self.death_rate.grid(row=1, column=1,pady=5,padx=15)
  
    self.albedo_black = Scale(self.container, from_=0, to=1,digits = 3, resolution = 0.01, orient=HORIZONTAL,label = 'Black Albedo',font=self.times18,width=40)
    self.albedo_black.grid(row=1, column=0,pady=5,padx=15)

    self.albedo_white = Scale(self.container, from_=0, to=1,digits = 3, resolution = 0.01, orient=HORIZONTAL,label = 'White Albedo',font=self.times18,width=40)
    self.albedo_white.grid(row=2, column=0,pady=5,padx=15)

    self.q_value = Scale(self.container, from_=0, to=3,digits = 3, resolution = 0.1, orient=HORIZONTAL,label = 'q(constant)',font=self.times18,width=40)
    self.q_value.grid(row=2, column=1,pady=5,padx=15)

    self.max_L2 = Scale(self.container, from_=1.6, to=3,digits = 4, resolution=0.01, orient=HORIZONTAL,label='Luminosity',font=self.times18,width=40)
    self.max_L2.grid(row=3, column = 0,pady=5,padx=15)

    self.button5 = Button(self.container,text='Black',command=self.area_black_daisies,font=self.times18)
    self.button5.grid(row=3,column=1,pady=5,padx=15)

    self.button6 = Button(self.container,text='White',command=self.area_white_daisies,font=self.times18)
    self.button6.grid(row=4,column=1,pady=5,padx=15)

    self.button6 = Button(self.container,text='Both',command=self.both_daisies,font=self.times18)
    self.button6.grid(row=5,column=1,pady=5,padx=15)

    self.return_button = Button(self.container,text='Exit',command = self.exit,font=self.times18)
    self.return_button.grid(row=4,column=0,pady=15)

    
  def set_death_rate(self):
    self.g = self.death_rate.get()

  def set_albedo_black(self):
    self.Ab = self.albedo_black.get()

  def set_albedo_white(self):
    self.Aw = self.albedo_white.get()

  def set_q_value(self):
    self.q = self.q_value.get() * (10**-9)

  def set_max_L(self):
    self.max_L = self.max_L2.get()
    self.L_step = (self.max_L - 0.6) / 1000
    self.L = np.arange(0.6,self.max_L,self.L_step)

  def albedo_calculator(self):
    self.two_species_area()
    self.set_albedo_black()
    self.set_albedo_white()
    self.A_list = []
    b_list = []
    w_list = []
    for i in range(len(self.ab_2_list)):
      b = self.ab_2_list[i] * self.Ab
      b_list.append(b)
    for i in range(len(self.aw_2_list)):
      w = self.aw_2_list[i] * self.Aw
      w_list.append(b)
    for i in range(len(b_list)):
      A = ((1-b_list[i]-w_list[i])*0.5) + (b_list[i]*self.Ab) + (w_list[i]*self.Aw)
      self.A_list.append(A)     

  def area_black_daisies(self):
    self.set_albedo_black()
    self.albedo_calculator()
    self.set_death_rate()
    self.set_q_value()
    self.set_max_L()
    self.mean_temp_black()
    q_prime = self.q/4 * ((273+22.5)**3)
    self.Tb2_list = []
    for i in range(len(self.Te_b_list)):
      Tb2 = q_prime*(self.A_list[i] - self.Ab) + self.Te_b_list[i]
      self.Tb2_list.append(Tb2)
    self.Tb_list = self.Tb2_list
    self.black_area()
    self.sorting_ab()
    self.plot(self.L[:len(self.ab_list)],self.ab_list,'Luminosity','Area')
    

  def area_white_daisies(self):
    self.set_albedo_white()
    self.albedo_calculator()
    self.set_death_rate()
    self.set_q_value()
    self.set_max_L()
    self.mean_temp_white()
    q_prime = self.q/4 * ((273+22.5)**3)
    self.Tw2_list = []
    for i in range(len(self.Te_w_list)):
      Tw2 = q_prime*(self.A_list[i] - self.Aw) + self.Te_w_list[i]
      self.Tw2_list.append(Tw2)
    self.Tw_list = self.Tw2_list
    self.white_area()
    self.sorting_aw()
    self.plot(self.L[:len(self.aw_list)],self.aw_list,'Luminosity','Area')

  def both_daisies(self):
   self.set_albedo_white()
   self.albedo_calculator()
   self.set_death_rate()
   self.set_q_value()
   self.set_max_L()
   self.mean_temp_black()
   q_prime = self.q/4 * ((273+22.5)**3)
   self.Tb2_list = []
   for i in range(len(self.Te_b_list)):
     Tb2 = q_prime*(self.A_list[i] - self.Ab) + self.Te_b_list[i]
     self.Tb2_list.append(Tb2)
   self.Tb_list = self.Tb2_list
   self.black_area()
   self.sorting_ab()
   aw_list = []
   for i in range(len(self.ab_list)):
     aw = 0.7 - self.ab_list[i]
     aw_list.append(aw)
   self.mul_plots(self.L[:len(self.ab_list)],self.ab_list,aw_list,'total area','luminosity')

  def exit(self):
      self.clearFrame()
      GUI(root) 
    

daisy = main()
root = Tk()
daist_GUI = GUI(root)

