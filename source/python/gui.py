import customtkinter
import api

small_font = ('Arial', 16)
big_font = ('Arial', 22)

class ResultWindow(customtkinter.CTkToplevel):
    def __init__(self, results, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Wyniki Obliczeń - Metoda Bairstowa")
        # Zwiększamy rozmiar okna, aby pomieściło długie liczby
        self.geometry("700x400")
        self.minsize(650, 300)

        # Nagłówek okna
        self.title_label = customtkinter.CTkLabel(
            self, 
            text="Znalezione Pierwiastki", 
            font=customtkinter.CTkFont(size=34, weight="bold")
        )
        self.title_label.pack(pady=(20, 10))

        # Dodajemy scrollowaną ramkę na wypadek wielu pierwiastków
        self.scrollable_frame = customtkinter.CTkScrollableFrame(self, fg_color="transparent")
        self.scrollable_frame.pack(padx=20, pady=(0, 20), fill="both", expand=True)

        # Wyświetlanie sformatowanych danych
        self.display_results(results)

    def display_results(self, results_text):
        # 1. Definiujemy czcionki na samym początku (używamy niezawodnych krotek)
        mono_font = ("Consolas", 20)
        bold_font = ("Arial", 30, "bold")

        lines = results_text.replace('?', '\n').split('\n')
        roots = [line.strip() for line in lines if "Root:" in line]

        # 2. Naprawiamy fallback - dodajemy dużą czcionkę i dopisek DEBUG!
        if not roots:
            fallback_label = customtkinter.CTkLabel(
                self.scrollable_frame, 
                text=f"{results_text}", 
                wraplength=600,
                font=mono_font  # Teraz awaryjny tekst TEŻ będzie wielki
            )
            fallback_label.pack(padx=10, pady=10)
            return

        # 3. Główna pętla
        for root_str in roots:
            card_frame = customtkinter.CTkFrame(self.scrollable_frame, corner_radius=8)
            card_frame.pack(padx=10, pady=8, fill="x")
            
            parts = root_str.split(" ", 2) 
            
            if len(parts) >= 3:
                header_text = f"{parts[0]} {parts[1]}"  # "Root: X"
                value_text = parts[2]                   # Wartości w nawiasach
                
                header_label = customtkinter.CTkLabel(card_frame, text=header_text, font=bold_font)
                header_label.pack(anchor="w", padx=15, pady=(10, 0))
                
                value_label = customtkinter.CTkLabel(card_frame, text=value_text, font=mono_font, justify="left")
                value_label.pack(anchor="w", padx=15, pady=(2, 10))
            else:
                single_label = customtkinter.CTkLabel(card_frame, text=root_str, font=mono_font)
                single_label.pack(padx=15, pady=10, anchor="w")

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
        print(self.counting_frame.get_value())
        print(self.entry.get())

        self.api.send_command(self.counting_frame.get_value())
        self.api.send_command(self.entry.get())
        
        # Pobieranie danych 
        # (Uwaga: jeśli po dodaniu 'cout << endl;' w C++ wysyłana jest tylko jedna linia, 
        # drugie wywołanie get_data() zawiesi program. Zostawiłem oba, ale miej to na uwadze).
        self.entry_text = self.api.get_data()
        self.entry_text = self.api.get_data()
        
        # Uproszczone wypisywanie do konsoli (zastępuje pętlę for i instrukcje if)
        self.entry_text = self.entry_text.replace('?', '\n')
        if self.entry_text:
            print(self.entry_text)

        # Logika wywołania nowego okienka z wynikami
        if self.toplevel_window is None or not self.toplevel_window.winfo_exists():
            self.toplevel_window = ResultWindow(results=self.entry_text)
        else:
            # Okienko istnieje - zamykamy je i tworzymy nowe, aby załadować nowe wyniki
            self.toplevel_window.destroy()
            self.toplevel_window = ResultWindow(results=self.entry_text)
            
        # Przeniesienie focusu na nowo otwarte okno
        self.toplevel_window.focus()

        
        


def main():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    main()