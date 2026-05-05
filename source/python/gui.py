import customtkinter
import api

small_font = ('Arial', 16)
big_font = ('Arial', 22)

class ResultWindow(customtkinter.CTkToplevel):
    def __init__(self, results, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Wyniki Obliczeń")
        self.geometry("300x200")

        # Wyświetlanie przekazanych danych
        self.label = customtkinter.CTkLabel(self, text=f"Wynik operacji: {results}")
        self.label.pack(padx=20, pady=20)


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

        self.api = api.API()

        self.toplevel_window = None

        self.title("my app")
        self.geometry("1050x600")
        self.grid_columnconfigure((0, 1, 2), weight=1)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        self.entry_text = "Wielomian w formacie: a_n,n a_(n-1),n-1 ... a_0,0 lub [a_n,b_n],n [a_(n-1),b_(n-1)],n-1 ... [a_0,b_0],0"

        self.title_text = customtkinter.CTkLabel(self, 
                                                 text="Kalkulator pierwiastków zespolonych wielomianu (metoda Bairstowa)",
                                                 fg_color="gray30", corner_radius=6,
                                                 font=big_font,height=40)
        self.title_text.grid(row=0, column=0, columnspan=3, padx=20, pady=20, sticky="nwes")

        self.counting_frame = RadiobtnFrame(self,
            "Typ zliczania", 
            [
                "zwykłe (arytmetyka zmiennoprzecinkowa)", 
                "przedziałowe dla danych rzeczywistych (arytmetyka przedziałowa)",
                "przedziałowe dla danych przedziałowych (arytmetyka przedziałowa)"
            ])
        
        self.counting_frame.grid(row=1, column=1, columnspan=2, padx=20, pady=20, sticky="nw")

        self.input_frame = customtkinter.CTkFrame(self)
        self.input_frame.grid(row=1, column=0, padx=20, pady=20, sticky="nw")

        self.run_button = customtkinter.CTkButton(self.input_frame, text="Uruchom", command=self.button_callback,font=big_font,height=80)
        self.run_button.grid(row=0, column=0, padx=40, pady=40, sticky="nwes")

        self.entry = customtkinter.CTkEntry(self, placeholder_text=self.entry_text,height=60,font=big_font)
        self.entry.grid(row=2, column=0,columnspan=3 ,padx=20, pady=20, sticky="new")

    def button_callback(self):
        print("Counting type:", self.counting_frame.get_value())
        print("Polynomial:", self.entry.get())

        self.api.send_command(self.counting_frame.get_value())
        #self.api.send_command(self.input_type_frame.get_value())
        self.api.send_command(self.entry.get())

        self.entry_text = self.api.get_data()
        print(self.entry_text)

        if self.toplevel_window is None or not self.toplevel_window.winfo_exists():
            self.toplevel_window = ResultWindow(results=self.entry_text)
        else:
            self.toplevel_window.focus()

        
        


def main():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    main()