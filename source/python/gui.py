import customtkinter

small_font = ('Arial', 16)
big_font = ('Arial', 22)

class RadiobtnFrame(customtkinter.CTkFrame):
    def __init__(self, master, title, chechbox_data : list[str]):
        super().__init__(master)
        self.variable = customtkinter.StringVar(value="")

        self.title = customtkinter.CTkLabel(self, text=title, fg_color="gray30", corner_radius=6,font=big_font)
        self.title.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="nwes")

        self.checkboxes = []

        for i in range(len(chechbox_data)):
            self.checkboxes.append(customtkinter.CTkRadioButton(self, text=chechbox_data[i],
                                                                variable=self.variable, value=i,font=small_font))  
            self.checkboxes[i].grid(row=i+1, column=0, padx=10, pady=(10, 10), sticky="nwes")
    
    def get_value(self):
        return self.variable.get()

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title("my app")
        self.geometry("1050x600")
        self.grid_columnconfigure((0, 1, 2), weight=1)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        self.title_text = customtkinter.CTkLabel(self, 
                                                 text="Kalkulator pierwiastków zespolonych wielomianu (metoda Bairstowa)",
                                                 fg_color="gray30", corner_radius=6,
                                                 font=big_font,height=40)
        self.title_text.grid(row=0, column=0, columnspan=3, padx=20, pady=20, sticky="nwes")

        self.arithemetic_frame = RadiobtnFrame(self,
            "Typ arytmetyki", 
            [
                "zwykła zmiennopozycyjna", 
                "przedziałowa zmiennopozycyjna"
            ])
        
        self.arithemetic_frame.grid(row=1, column=1, padx=20, pady=20, sticky="nw")

        self.counting_frame = RadiobtnFrame(self,
            "Typ zliczania", 
            [
                "zwykłe", 
                "przedziałowe dla danych rzeczywistych",
                "przedziałowe dla danych przedziałowych"
            ])
        
        self.counting_frame.grid(row=1, column=2, padx=20, pady=20, sticky="nw")

        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=20, pady=20, sticky="nw")

        self.run_button = customtkinter.CTkButton(self.input_frame, text="Uruchom", command=self.button_callback,font=big_font,height=80)
        self.run_button.grid(row=0, column=0, padx=40, pady=40, sticky="nwes")

        self.input_type_frame = RadiobtnFrame(self.input_frame,
            "Typ danych wejściowych", 
            [
                "tekst", 
                "czytanie z pliku"
            ])
        
        self.input_type_frame.grid(row=1, column=0, padx=20, pady=20, sticky="nw")


        self.entry = customtkinter.CTkEntry(self, placeholder_text="Wprowadź wielomian",height=60,font=big_font)
        self.entry.grid(row=2, column=0,columnspan=3 ,padx=20, pady=20, sticky="new")

    def button_callback(self):
        print("Arithmetic type:", self.arithemetic_frame.get_value())
        print("Counting type:", self.counting_frame.get_value())
        print("Input type:", self.input_type_frame.get_value())
        print("Polynomial:", self.entry.get())


def main():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    main()