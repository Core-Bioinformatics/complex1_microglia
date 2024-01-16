import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('../processed/complex2_clusters.csv')
df = pd.melt(df, value_vars=['A-EAE_FALSE', 'A-EAE_TRUE', 'C-EAE_FALSE', 'C-EAE_TRUE',
       'Ctrl_FALSE', 'Ctrl_TRUE'], id_vars='Unnamed: 0')

df['variable'] = df['variable'].apply(lambda x: '{}\n({})'.format(x.split('_')[0], 'Not in C4' if x.endswith('FALSE') else 'In C4'))

sns.barplot(data=df, x='variable', y='value', palette='Set2', errwidth=0, capsize=0.3)
plt.title('Average CII expression by state and cluster')
plt.ylabel('Avg. expression')
plt.xlabel('Class')
plt.grid(alpha=0.2)
