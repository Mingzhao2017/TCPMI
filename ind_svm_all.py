import sklearn
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, average_precision_score
from sklearn.svm import SVC
import pandas as pd

data_names = ['Lv']
all_results = []

for dataname in data_names:
    print(dataname)
    matrix_data = np.loadtxt('./data_coded/data_' + dataname + '_Tri_BiPSDP.csv', delimiter=' ')
    label_data = np.loadtxt('./data_coded/label_' + dataname + '.csv', delimiter=' ')

    # 找到正类和负类索引
    positive_indices = np.where(label_data == 1)[0]
    negative_indices = np.where(label_data != 1)[0]

    # 按规则选取 1/20 正类
    sampled_positive_indices = positive_indices[::3]  # 每隔5个取1条

    # 训练集索引（1/20正类 + 所有负类）
    train_indices = np.concatenate([sampled_positive_indices, negative_indices])

    # 测试集索引（剩余正类）
    test_indices = np.setdiff1d(positive_indices, sampled_positive_indices)

    # 划分数据
    train_matrix = matrix_data[train_indices]
    train_label = label_data[train_indices]
    test_matrix = matrix_data[test_indices]
    test_label = label_data[test_indices]

    print(train_matrix.shape, train_label.shape, test_matrix.shape, test_label.shape)

    # 初始化结果列表
    results = []

    # 训练 SVM
    svm_classifier = SVC(kernel='rbf', probability=True, C=0.5)
    svm_classifier.fit(train_matrix, train_label)

    # 预测测试集
    test_pred = svm_classifier.predict(test_matrix)

    # 计算指标
    accuracy = accuracy_score(test_label, test_pred)
    tn, fp, fn, tp = confusion_matrix(test_label, test_pred).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    mcc_value = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) if (tp + fp) * (
                tp + fn) * (tn + fp) * (tn + fn) > 0 else 0
    auroc_value = roc_auc_score(test_label, svm_classifier.predict_proba(test_matrix)[:, 1])
    auprc_value = average_precision_score(test_label, svm_classifier.predict_proba(test_matrix)[:, 1])

    # 添加结果
    results.append({
        'Fold': 10,
        'Accuracy': accuracy,
        'Sensitivity': sensitivity,
        'Specificity': specificity,
        'MCC': mcc_value,
        'AUROC': auroc_value,
        'AUPRC': auprc_value,
    })

    avg_results = {key: np.mean([result[key] for result in results]) for key in results[0] if key != 'Fold'}
    avg_results['Dataset'] = 'DNA_6mA_' + dataname + '_' + str(train_matrix.shape[0]) + ' ' + 'DNA_6mA_' + '_' + str(
        test_matrix.shape[0])
    all_results.append(avg_results)

    # ===== 新增：保存预测结果到 CSV =====
    pred_df = pd.DataFrame({
        'TrueLabel': test_label.astype(int),
        'PredLabel': test_pred.astype(int)
    })
    pred_df.to_csv('pred.csv', index=False)
    print("预测结果已保存到 pred.csv")

# 保存最终评估结果
all_results_df = pd.DataFrame(all_results)
all_results_df.to_excel('model_ind_results_svm_0.5_Lv_Lv_ind.xlsx', index=False)
