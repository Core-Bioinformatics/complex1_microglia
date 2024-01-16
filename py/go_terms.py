import pandas as pd
import seaborn as sns

df_go = {}

for cluster in range(3):
    df_go[cluster] = pd.read_csv('../processed/go/subclusters/go_{}.csv'.format(cluster))
    df_go[cluster] = df_go[cluster].drop_duplicates(subset='intersection').head(3)
    df_go[cluster]['negative_log_pval'] = -np.log(df_go[cluster]['p_value'])
    df_go[cluster]['name'] = df_go[cluster].apply(lambda x: '{}\n({})'.format(x['term_name'], x['term_id']), axis=1)

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, facecolor='white', alpha=0, figsize=(10, 5))

for row, (cluster, temp_df) in enumerate(df_go.items()):
    ax = axes[row]
    sns.barplot(data=temp_df[['name', 'negative_log_pval']], x='negative_log_pval', y='name', palette='Set2', ax=ax)
    ax.set_ylabel('')
    ax.set_title(cluster)
    if row != 2:
        ax.set_xlabel('')
    else:
        ax.set_xlabel('Negative log p-value')
    ax.grid(alpha=0.2)
        
plt.suptitle('GO terms for subclusters')
plt.tight_layout()
plt.savefig('../figs/goterms.pdf', bbox_inches="tight")
