from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import mean_squared_error, classification_report, ConfusionMatrixDisplay, confusion_matrix, f1_score
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
import pandas as pd


def clasificador(metadata, df):

    variables = ['ER', 'HER2', 'B-R grade', 'node status', 'tumor type']
    df_metrics = pd.DataFrame(columns=['VAR', 'max_depth', 'n_estimators', 'acc_train', 'acc_test', 'f1_score', 'rmse'])

    for var in variables:
        # Se obtiene histogrma de cada fenotipo
        plt.figure()
        plt.figure(figsize=(12, 6))
        categoria_counts = metadata[var].value_counts()
        categoria_counts.plot(kind='bar')
        plt.savefig("./results/clasificador/histograma_"+str(var)+".jpg")

        # Se binarizan las clases de cada variable objetivo
        objetivo = OneHotEncoder().fit_transform(metadata[[var]]).toarray()

        # Se crean los conjuntos de train y test
        X_train, X_test, y_train, y_test = train_test_split(df, objetivo, test_size = 0.3, random_state=5)

        # Se localizan los parámetros óptimos para el modelo
        ran_forest = RandomForestClassifier()
        param_grid = {'max_depth': range(6,13), 'n_estimators': [10, 50, 100, 200]}

        grid = GridSearchCV(ran_forest, param_grid=param_grid, cv=5)
        grid.fit(X_train, y_train)

        # Se crea el modelo con los parámetros obtenidos        
        rf_best = RandomForestClassifier(max_depth= grid.best_params_['max_depth'], n_estimators= grid.best_params_['n_estimators'], class_weight="balanced")
        rf_best.fit(X_train, y_train)

        # Se realiza la predicción y se evalua el modelo
        y_pred = rf_best.predict(X_test)
        rmse = mean_squared_error(y_true  = y_test.argmax(axis=1), y_pred  = y_pred.argmax(axis=1), squared = False)
        acc_train = rf_best.score(X_train, y_train)
        acc_test = rf_best.score(X_test, y_test)

        # Se guardan los resultados
        df_newrow = pd.DataFrame({'VAR': [var], 'max_depth': [grid.best_params_['max_depth']], 'n_estimators': [grid.best_params_['n_estimators']], 'acc_train': [acc_train], 'acc_test': [acc_test], 'rmse': [rmse]})
        df_metrics = pd.concat([df_metrics, df_newrow])
        
        # Report 
        report = classification_report(y_test.argmax(axis=1), y_pred.argmax(axis=1), digits=3, output_dict=True)
        df_report = pd.DataFrame(report).transpose()
        df_report.to_csv("./results/clasificador/report_"+str(var)+".csv")

        # Matriz de confusión
        cnf_matrix = confusion_matrix(y_test.argmax(axis=1), y_pred.argmax(axis=1))
        cm_label = sorted(metadata[var].unique().tolist())
        plt.figure()
        plt.figure(figsize=(12, 6))
        cm_plot = ConfusionMatrixDisplay(confusion_matrix = cnf_matrix, display_labels = cm_label)
        cm_plot.plot()
        plt.savefig("./results/clasificador/matriz_confusion_"+str(var)+".jpg")
    
    df_metrics.to_csv("./results/clasificador/resultados.csv")
        