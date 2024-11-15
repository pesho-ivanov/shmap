import marimo

__generated_with = "0.9.17"
app = marimo.App()


@app.cell
def __():
    from readpaf import parse_paf
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import marimo as mo
    return mo, np, parse_paf, pd, plt


@app.cell
def __(parse_paf):
    handle = open('evals/simple_test_gt.paf')
    df = parse_paf(handle, dataframe=True)
    return df, handle


@app.cell
def __(df, np, parse_paf, plt):
    def plot_corr(f):
        handle = open(f)
        df = parse_paf(handle, dataframe=True)

        plt.xlim(0.5, 0.9)
        plt.ylim(0.5, 0.9)
        plt.scatter(df.gt_J + np.random.normal(0, 0.01, size=len(df.gt_J)),
                    df.gt_C + np.random.normal(0, 0.01, size=len(df.gt_C)),
                    alpha=0.1)
        plt.xlabel('gt_J')
        plt.ylabel('gt_C')
        plt.show()

    # get a subset of columns of df
    df2 = df[['J', 'J2', 'gt_J', 'gt_C', 'gt_C_bucket']]
    return df2, plot_corr


@app.cell
def __(df2, pd, plt):
    axes = pd.plotting.scatter_matrix(df2, alpha=0.1)
    for i in range(len(axes)):
        for j in range(len(axes[i])):
            axes[i, j].set_xticks([0.55, 0.85])
            axes[i, j].set_yticks([0.55, 0.85])
            if i != j:  # Skip diagonal plots
                axes[i, j].set_xlim(0.5, 0.9)
                axes[i, j].set_ylim(0.5, 0.9)
    plt.gca()
    return axes, i, j


@app.cell
def __(plot_corr):
    plot_corr('chm13_gt.paf')
    return


if __name__ == "__main__":
    app.run()
