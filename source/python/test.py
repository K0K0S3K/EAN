import customtkinter as ctk

class ResultWindow(ctk.CTkToplevel):
    def __init__(self, results, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Wyniki Obliczeń")
        self.geometry("300x200")

        # Wyświetlanie przekazanych danych
        self.label = ctk.CTkLabel(self, text=f"Wynik operacji: {results}")
        self.label.pack(padx=20, pady=20)

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Panel Sterowania")
        self.geometry("400x300")

        self.entry = ctk.CTkEntry(self, placeholder_text="Wpisz dane...")
        self.entry.pack(pady=20)

        self.button = ctk.CTkButton(self, text="Oblicz i pokaż wyniki", command=self.run_engine)
        self.button.pack(pady=10)
        
        self.toplevel_window = None

    def run_engine(self):
        data = self.entry.get()
        
        # Symulacja silnika obliczeniowego
        processed_data = f"Przetworzono: {data.upper()}"
        
        # Sprawdzanie, czy okno już istnieje, żeby nie otwierać kilku naraz
        if self.toplevel_window is None or not self.toplevel_window.winfo_exists():
            self.toplevel_window = ResultWindow(results=processed_data)
        else:
            self.toplevel_window.focus()  # Jeśli istnieje, wysuń na wierzch

if __name__ == "__main__":
    app = App()
    app.mainloop()