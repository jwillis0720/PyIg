from Tkinter import *
import tkFileDialog 

class pyigblast_gui():
	def __init__(self, parent):

		

		#------ constants for controlling layout of buttons ------
		button_width = 6
		button_padx = "2m"
		button_pady = "1m"
		buttons_frame_padx =  "3m"
		buttons_frame_pady =  "2m"
		buttons_frame_ipadx = "3m"
		buttons_frame_ipady = "1m"
		# -------------- end constants ----------------


		# set up Tkinter variables, to be controlled by the radio buttons
		self.button_name = StringVar()
		self.button_name.set("C")


		# -------------- end constants ----------------

		self.myParent = parent
		self.myParent.geometry("750x750")

		### Our topmost frame is called myContainer1
		self.container = Frame(parent) ### Frame is a class in tkinter
		self.container.pack(expand=YES, fill=BOTH)

		# control frame - all the options for pyigblast
		self.control_frame = Frame(self.container) ###
		self.control_frame.pack(side=TOP, expand=YES,  padx=10, pady=5, ipadx=5, ipady=5)

		# inside control_frame we create a header label
		# and a buttons_frame at the top,
		# and demo_frame at the bottom

		Label(self.control_frame, text="PyIgBLAST", justify="center",font=("Arial", 24)).pack(side=TOP, anchor=NW)
		Label(self.control_frame, text="     Jordan Willis", justify="center",font=("Arial", 16)).pack(side=TOP, anchor=NW)
		
		# buttons frame
	
		self.buttons_frame = Frame(self.control_frame) ###
		
		self.buttons_frame.pack(side=LEFT, expand=NO, fill=BOTH, ipadx=5, ipady=5)
		self.fasta_file = StringVar()
		self.fasta_file_button = Button(self.buttons_frame,text="FASTA file",command=self.askopenfile).pack()
		#self.fasta_file = self.fasta_file_button.get()





		# output frame
		#self.demo_frame = Frame(self.container) ###
		
		#don't pack anything in the frame until the output 
		#self.demo_frame.pack(side=BOTTOM, expand=YES, fill=BOTH)


		self.cancelButtonFrame = Frame(self.container)
		self.cancelButtonFrame.pack(side=BOTTOM, expand=YES, anchor=SW)
		self.cancelButton = Button(self.cancelButtonFrame,
			text="Cancel", background="red",
			width=button_width,
			padx=button_padx,
			pady=button_pady
			)
		
		self.cancelButton.pack(side=BOTTOM, anchor=S)
		self.cancelButton.bind("<Button-1>", self.cancelButtonClick)
		self.cancelButton.bind("<Return>", self.cancelButtonClick)
		
	def askopenfile(self):
		return tkFileDialog.askopenfilename(defaultextension='.fasta',initialdir='.')

	def cancelButtonClick(self, event):
		self.myParent.destroy()

		### Inside the demo frame, we create top_frame and bottom_frame.
		### These will be our demonstration frames.
		# top frame
		# self.top_frame = Frame(self.demo_frame)
		# self.top_frame.pack(side=TOP, expand=YES, fill=BOTH)  ###

		# # bottom frame
		# self.bottom_frame = Frame(self.demo_frame,
		#     borderwidth=5, 	relief=RIDGE,
		#     height=50,
		#     bg="cyan",
		#     ) ###
		# self.bottom_frame.pack(side=TOP, fill=X)


		### Now we will put two more frames, left_frame and right_frame,
		### inside top_frame.  We will use HORIZONTAL (left/right)
		### orientation within top_frame.

		# left_frame
		# self.left_frame = Frame(self.top_frame,	background="red",
		#     borderwidth=5, 	relief=RIDGE,
		#     width=50,
		#     ) ###
		# self.left_frame.pack(side=LEFT, expand=NO, fill=Y)


		# ### right_frame
		# self.right_frame = Frame(self.top_frame, background="tan",
		#     borderwidth=5, 	relief=RIDGE,
		#     width=250
		#     )
		# self.right_frame.pack(side=RIGHT, expand=YES, fill=BOTH)


		# now put a button in each of the interesting frames
		# button_names = ["A", "B", "C"]
		# side_options = [LEFT, TOP, RIGHT, BOTTOM]
		# fill_options = [X, Y, BOTH, NONE]
		# expand_options = [YES, NO]
		# anchor_options = [NW, N, NE, E, SE, S, SW, W, CENTER]


		# self.buttonA = Button(self.bottom_frame, text="A")
		# self.buttonA.pack()
		# self.buttonB = Button(self.left_frame, text="B")
		# self.buttonB.pack()
		# self.buttonC = Button(self.right_frame, text="C")
		# self.buttonC.pack()
		# self.button_with_name = {"A":self.buttonA, "B":self.buttonB, "C":self.buttonC}

		# now we some subframes to the buttons_frame
		# self.button_names_frame   = Frame(self.buttons_frame, borderwidth=5)
		# self.side_options_frame   = Frame(self.buttons_frame, borderwidth=5)
		# self.fill_options_frame   = Frame(self.buttons_frame, borderwidth=5)
		# self.expand_options_frame = Frame(self.buttons_frame, borderwidth=5)
		# self.anchor_options_frame = Frame(self.buttons_frame, borderwidth=5)

		# self.button_names_frame.pack(  side=LEFT, expand=YES, fill=Y, anchor=N)
		# self.side_options_frame.pack(  side=LEFT, expand=YES, anchor=N)
		# self.fill_options_frame.pack(  side=LEFT, expand=YES, anchor=N)
		# self.expand_options_frame.pack(side=LEFT, expand=YES, anchor=N)
		# self.anchor_options_frame.pack(side=LEFT, expand=YES, anchor=N)

		# Label(self.button_names_frame, text="\nButton").pack()
		# Label(self.side_options_frame, text="Side\nOption").pack()
		# Label(self.fill_options_frame, text="Fill\nOption").pack()
		# Label(self.expand_options_frame, text="Expand\nOption").pack()
		# Label(self.anchor_options_frame, text="Anchor\nOption").pack()

		for option in button_names:
			button = Radiobutton(self.button_names_frame, text=str(option), indicatoron=1,
				value=option, command=self.button_refresh, variable=self.button_name)
			button["width"] = button_width
			button.pack(side=TOP)

		for option in side_options:
			button = Radiobutton(self.side_options_frame, text=str(option), indicatoron=0,
				value=option, command=self.demo_update, variable=self.side_option)
			button["width"] = button_width
			button.pack(side=TOP)

		for option in fill_options:
			button = Radiobutton(self.fill_options_frame, text=str(option), indicatoron=0,
				value=option, command=self.demo_update, variable=self.fill_option)
			button["width"] = button_width
			button.pack(side=TOP)

		for option in expand_options:
			button = Radiobutton(self.expand_options_frame, text=str(option), indicatoron=0,
				value=option, command=self.demo_update, variable=self.expand_option)
			button["width"] = button_width
			button.pack(side=TOP)

		for option in anchor_options:
			button = Radiobutton(self.anchor_options_frame, text=str(option), indicatoron=0,
				value=option, command=self.demo_update, variable=self.anchor_option)
			button["width"] = button_width
			button.pack(side=TOP)




		self.cancelButton.bind("<Button-1>", self.cancelButtonClick)
		self.cancelButton.bind("<Return>", self.cancelButtonClick)

		# set up the buttons in their initial position
		self.demo_update()


	def button_refresh(self):
		button = self.button_with_name[self.button_name.get()]
		properties = button.pack_info()
		self.fill_option.set  (  properties["fill"] )
		self.side_option.set  (  properties["side"] )
		self.expand_option.set(  properties["expand"] )
		self.anchor_option.set(  properties["anchor"] )


	def demo_update(self):
		button = self.button_with_name[self.button_name.get()]
		button.pack(fill=self.fill_option.get()
			, side=self.side_option.get()
			, expand=self.expand_option.get()
			, anchor=self.anchor_option.get()
			)




root = Tk() #Root houses the main frame
myapp = pyigblast_gui(root)
root.mainloop()