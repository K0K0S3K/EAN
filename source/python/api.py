import subprocess
import os
import sys

class API:
    def __init__(self):
        self.algorithm_path = "./source/cpp/main"
        self.process = subprocess.Popen(
                self.algorithm_path,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
                )
        #self.run_algorithm()

        print(self.get_data())


    def send_command(self, cmd):
        self.process.stdin.write(f"{cmd}\n")
        self.process.stdin.flush()

    def get_data(self):
        return self.process.stdout.readline().strip()

    def run_algorithm(self):
        try:
            if self.process.returncode == 0:
                print("Algorithm output:")
                print(self.process.stdout)
            else:
                print("Algorithm error:")
                print(self.process.stderr)
        except Exception as e:
            print(f"An error occurred while running the algorithm: {e}")