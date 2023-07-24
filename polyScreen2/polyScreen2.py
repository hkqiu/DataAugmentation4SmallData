#coding=utf-8

import os,sys
import tkinter as tk
import webbrowser

dir_path = os.path.dirname(__file__)
dependency_path = os.path.join(dir_path, 'dependency')
sys.path.insert(0, dependency_path)
import deepchem as dc
import pandas as pd
from rdkit import Chem


class PolyScreen:
    def __init__(self, master, width=500, height=400):
        self.master = master  # master
        self.master.title("polyScreen")
        self.w = width
        self.h = height
        self.__set_menu()
        self.path = tk.StringVar()  # For Properties Prediction
        self.property = tk.StringVar()  # For Inverse Design
        self.dialogue = tk.StringVar()  # For PolyChat

        self.job_path = os.path.dirname(os.path.dirname(__file__))  # path of current directory
        self.MODEL_DIR =  os.path.join(self.job_path, 'model/ToObtainPredictedValuesofPolyimides')

        self.pi_alternative_path = 'PATH/TO/ALTERNATIVE/csv'

        self.button = tk.Button(self.master, text="Properties Prediction", font=('Times New Roman', 25), command=self.__properties_prediction)
        self.button.pack(pady=10)
        self.button = tk.Button(self.master, text="Inverse Design", font=('Times New Roman', 25), command=self.__inverse_design)
        self.button.pack(pady=10)
        self.button = tk.Button(self.master, text="PolyChat", font=('Times New Roman', 25), command=self.__poly_chat)
        self.button.pack(pady=10)

        # author = """
        # ©2023 Sun's AI4P Workshop. All rights reserved.
        # Email: zysun@ciac.ac.cn
        # """

        text = tk.Text(root, height=3)
        text.pack(side=tk.BOTTOM)
        text.insert(tk.END, "©2023 Sun's AI4P Workshop. All rights reserved.\nEmail: zysun@ciac.ac.cn")

    def __center(self):
        ws = self.master.winfo_screenwidth()
        hs = self.master.winfo_screenheight()
        x = int((ws / 2) - (self.w / 2))
        y = int((hs / 2) - (self.h / 2))
        self.master.geometry('{}x{}+{}+{}'.format(self.w, self.h, x, y))

    def __set_menu(self):
        # Menu菜单
        menu = tk.Menu(self.master)
        self.master.config(menu=menu)  # 在窗口添加菜单栏
        combobox = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Links', menu=combobox)  # 在菜单栏添加下拉菜单

        combobox.add_command(label='Github',
                             command=lambda: webbrowser.open('https://github.com/HKQiu/DataEnhancement4SmallData'))
        combobox.add_command(label='Tutorial',
                             command=lambda: webbrowser.open('https://github.com/HKQiu/DataEnhancement4SmallData'))

    def __run(self):
        if self.path.get() != '':
            print("Filepath:", self.path.get())
        elif self.property.get() != '':
            print("Property:", self.property.get())
        elif self.dialogue.get() != '':
            print("Chat:", self.dialogue.get())

    def __properties_prediction(self):
        # 创建新窗口
        new_window = tk.Toplevel(self.master)
        new_window.geometry("800x480")
        new_window.title("Properties Prediction")

        lab1 = tk.Label(new_window, text='Path to the file:', font=('Times New Roman', 28))
        entry = tk.Entry(new_window, textvariable=self.path,
                         bd=10, font=('Times New Roman', 20), highlightcolor='Fuchsia', width=30)
        play = tk.Button(new_window, text="Run",
                         font=('Times New Roman', 20), fg='Purple', width=2, height=1, command=self.__run)

        lab1.grid(row=0, column=0)
        entry.grid(row=0, column=1)
        play.grid(row=0, column=3, ipadx=10, ipady=10)

    def __inverse_design(self):
        # 创建新窗口
        new_window = tk.Toplevel(self.master)
        new_window.geometry("800x480")
        new_window.title("Inverse Design")

        lab1 = tk.Label(new_window, text='PI with Tg over ', font=('Times New Roman', 28))
        entry = tk.Entry(new_window, textvariable=self.property,
                         bd=10, font=('Times New Roman', 20), highlightcolor='Fuchsia', width=5)
        lab2 = tk.Label(new_window, text='°C', font=('Times New Roman', 28))
        play = tk.Button(new_window, text="Run",
                         font=('Times New Roman', 20), fg='Purple', width=2, height=1, command=self.__run)

        lab1.grid(row=0, column=1)
        entry.grid(row=0, column=2)
        lab2.grid(row=0, column=3)
        play.grid(row=0, column=4, ipadx=10, ipady=10)

    def __poly_chat(self):
        # 创建新窗口
        new_window = tk.Toplevel(self.master)
        new_window.geometry("800x480")
        new_window.title("PolyChat")

        lab1 = tk.Label(new_window, text='Tell me what to do:', font=('Times New Roman', 28))
        entry = tk.Entry(new_window, textvariable=self.dialogue,
                         bd=10, font=('Times New Roman', 20), highlightcolor='Fuchsia', width=30)

        play = tk.Button(new_window, text="Run",
                         font=('Times New Roman', 20), fg='Purple', width=2,
                         height=1, command=self.__run)

        lab1.grid(row=0, column=0)
        entry.grid(row=0, column=1)
        play.grid(row=0, column=3, ipadx=10, ipady=10)

    def predict(self):
        DATASET_FILE = self.path.get()
        MODEL_DIR = self.MODEL_DIR

        featurizer = dc.feat.ConvMolFeaturizer()
        loader = dc.data.CSVLoader(tasks=[], feature_field="Smiles", featurizer=featurizer)
        testset = loader.create_dataset(DATASET_FILE, shard_size=10000)

        model = dc.models.GraphConvModel(1, mode="regression", model_dir=MODEL_DIR)
        model.restore()

        test_pred = model.predict(testset)
        pred = pd.DataFrame(test_pred)
        pred.columns = ['Pred/°C']
        smiles_df = pd.read_csv(DATASET_FILE)
        res = pd.concat([smiles_df, pred], axis=1)
        save_path = os.path.join(self.job_path, 'Temp', 'predict.csv')
        res.to_csv(save_path)

    def design(self, property_needed=None):
        pi_alternative_path = self.pi_alternative_path
        pi = pd.read_csv(pi_alternative_path)
        if not property_needed:
            property_needed = int(self.property.get())
            screend_pi = pi[pi['Pred'] >= property_needed]
            root_path = pi_alternative_path.split('.', 1)[0]
            save_path = os.path.join(self.job_path, 'Temp', 'design.csv')
            screend_pi.to_csv(save_path)
        else:
            property_needed = property_needed
            screend_pi = pi[pi['Pred'] >= property_needed]
            root_path = pi_alternative_path.split('.', 1)[0]
            save_path = os.path.join(self.job_path, 'Temp', 'chat_design.csv')
            screend_pi.to_csv(save_path)

    def chat(self):
        _chat = self.dialogue.get()
        token_list = _chat.split()

        inverse_design = ['design', 'screen']
        tg_predict = ['Tg', 'tg', 'Glass Transition temperature']

        for token in token_list:
            if token in inverse_design:
                task = 'design'
                print("Executing designing...")
                break
            elif token in tg_predict:
                pass  # 执行预测
                task = 'predict'
                break
            else:
                sun_poem = """
                水调歌头
                        ——送给热爱科研的同学们
                稳坐桌台前，十指飞一般。
                思如洪波细流，涌入脑海间。
                意欲弃之归去，心又有所不甘，感慨万万千。
                倾身若续茶，复转头钻研。
                写代码，编脚本，阅文献。
                时光飞逝，离别之际已在前。
                不负多年心血，能否终获答案？
                此时竟无言。
                归来常忆起，曾经苦与甜。
                """
                print(sun_poem)


        if task == 'predict':
            smiles = [i for i in token_list if self.__is_smiles(i)]
            smiles_df = pd.DataFrame(smiles, columns=['Smiles'])
            smiles_tem_dir = self.job_path
            if not os.path.exists(os.path.join(smiles_tem_dir, 'Temp')):
                os.mkdir(os.path.join(smiles_tem_dir, 'Temp'))
            smiles_df.to_csv(os.path.join(smiles_tem_dir, 'Temp', 'temp.csv'))

            featurizer = dc.feat.ConvMolFeaturizer()
            loader = dc.data.CSVLoader(tasks=[], feature_field="Smiles", featurizer=featurizer)
            testset = loader.create_dataset(os.path.join(smiles_tem_dir, 'Temp', 'temp.csv'),
                                            shard_size=10000)

            model = dc.models.GraphConvModel(1, mode="regression", model_dir=self.MODEL_DIR)
            model.restore()
            test_pred = model.predict(testset).tolist()
            res = {'Smiles': smiles, 'Pred/°C': test_pred}
            res = pd.DataFrame(res)
            res.to_csv(os.path.join(smiles_tem_dir, 'Temp', 'chat_predict.csv'))

        if task == 'design':
            property_needed = [int(i) for i in token_list if self.__transform_num(i)][0]
            self.design(property_needed)

    def __transform_num(self, string):
        """
        Determines whether a string can be converted to a number, returning True or False.
        """
        try:
            int(string)
            return True
        except ValueError:
            return False

    def __is_smiles(self, string):
        mol = Chem.MolFromSmiles(string)
        try:
            Chem.SanitizeMol(mol)
            return True
        except:
            return False

    def loop(self):
        # 禁止修改窗口大小
        self.master.resizable(False, False)
        # 窗口居中
        self.__center()
        self.master.mainloop()
        pass



root = tk.Tk()
app = PolyScreen(root)
app.loop()
if app.path.get() != '':
    app.predict()
elif app.property.get() != '':
    app.design()
elif app.dialogue.get() != '':
    app.chat()


if __name__ == '__main__':
    root = tk.Tk()
    app = PolyScreen(root)
    app.loop()
    if app.path.get() != '':
        app.predict()
    elif app.property.get() != '':
        app.design()
    elif app.dialogue.get() != '':
        app.chat()
