import customtkinter
import api

small_font = ('Arial', 16)
big_font = ('Arial', 22)

FONT_NORMAL = ("Segoe UI", 14)
FONT_BOLD = ("Segoe UI", 16, "bold")
FONT_MONO = ("Consolas", 15)

class ResultWindow(customtkinter.CTkToplevel):
    def __init__(self, results, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Wyniki Obliczeń - Metoda Bairstowa")
        self.geometry("1050x500")
        self.minsize(600, 400)
        
        self.after(100, self.lift)
        self.attributes("-topmost", True)

        self.header_label = customtkinter.CTkLabel(
            self, 
            text="Znalezione Pierwiastki", 
            font=customtkinter.CTkFont(size=28, weight="bold")
        )
        self.header_label.pack(pady=(20, 15))

        self.scrollable_frame = customtkinter.CTkScrollableFrame(
            self, 
            fg_color="transparent",
            label_text="Lista wyników"
        )
        self.scrollable_frame.pack(padx=20, pady=(0, 20), fill="both", expand=True)

        self.display_results(results)

    def display_results(self, results_text):
        clean_text = results_text.replace('?', '\n').strip()
        lines = [line.strip() for line in clean_text.split('\n') if line.strip()]

        if not lines:
            self._show_empty_state("Brak danych do wyświetlenia lub błąd formatu.")
            return

        roots_found = False
        for line in lines:
            if "root" in line.lower():
                self._create_root_card(line)
                roots_found = True
        
        if not roots_found:
            self._show_empty_state(clean_text)

    def _create_root_card(self, text):
        card = customtkinter.CTkFrame(self.scrollable_frame, corner_radius=10)
        card.pack(padx=5, pady=5, fill="x")

        if ":" in text:
            label_part, value_part = text.split(":", 1)
            
            lbl = customtkinter.CTkLabel(card, text=label_part.strip() + ":", font=FONT_BOLD, text_color="#1f6aa5")
            lbl.pack(anchor="w", padx=15, pady=(10, 0))
            
            val = customtkinter.CTkLabel(
                card, 
                text=value_part.strip(), 
                font=FONT_MONO, 
                wraplength=950, 
                justify="left"
            )
            val.pack(anchor="w", padx=15, pady=(0, 10))
        else:
            lbl = customtkinter.CTkLabel(card, text=text, font=FONT_NORMAL)
            lbl.pack(padx=15, pady=10, anchor="w")

    def _show_empty_state(self, message):
        error_label = customtkinter.CTkLabel(
            self.scrollable_frame, 
            text=message, 
            font=FONT_NORMAL,
            wraplength=500,
            text_color="gray"
        )
        error_label.pack(pady=40)

class RadiobtnFrame(customtkinter.CTkFrame):
    def __init__(self, master, title, chechbox_data : list[str]):
        super().__init__(master)
        self.variable = customtkinter.StringVar(value="0")

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
        # 1. Pobranie i wysłanie danych
        mode = self.counting_frame.get_value()
        poly_input = self.entry.get()

        if not poly_input.strip():
            print("Błąd: Pole wejściowe jest puste")
            return

        self.api.send_command(mode)
        self.api.send_command(poly_input)
    
        final_result = self.api.get_data()

        if self.toplevel_window is not None:
            self.toplevel_window.destroy()
            
        self.toplevel_window = ResultWindow(results=final_result)
        self.toplevel_window.focus()

        
        


def main():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    main()