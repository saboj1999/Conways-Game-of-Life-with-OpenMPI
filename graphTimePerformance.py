import matplotlib.pyplot as plt

def readFile(filename):
    _500x500_1000 = {}

    with open(filename) as file:
        print("Reading data from "+filename+" : \n")
        for line in file:
            time = float(line.strip().split(" ")[-1].strip())
            size = int(line.strip().split(" ")[-4].strip())
            xsize = int(line.strip().split(" ")[-5].strip())
            gen = int(line.strip().split(" ")[-3].strip())
            Processes = int(line.strip().split(" ")[6].strip())
            print(
                 "Time: "+str(time)+"\n"
                +"Size: "+str(size)+"x"+str(xsize)+"\n"
                +"Generations: "+str(gen)+"\n"
                +"Processes: "+str(Processes)+"\n"
            )
            _500x500_1000[Processes] = time      

    data = {"1024x1024 - 1000 Generations": _500x500_1000}     

    return data

def plotEfficiency(data):
    for key, value in data.items():
        plt.plot(list(value.keys()), list(value.values()), color='darkblue')
        plt.xscale('log')
        plt.xticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
        plt.xlabel('Processes')
        plt.ylabel('Time (seconds)')
        plt.title(key)
        plt.show() 

def GenerateSpeedUpData(data):
    # data is dict of Processes:time
    baseTime = data[1]
    for key, value in data.items():
        data[key] = baseTime / value
    
    return data

def plotSpeedUp(data):
    speedUpData = GenerateSpeedUpData(data["1024x1024 - 1000 Generations"])
    plt.plot(list(speedUpData.keys()), list(speedUpData.values()), color='darkred')
    plt.xscale('log')
    plt.xticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
    plt.xlabel('Processes')
    plt.ylabel('Speed-Up')
    plt.title("1024x1024 - 1000 Generations")
    plt.show() 

def main():
    data = readFile("time_average_outputs.txt")
    plotEfficiency(data)
    plotSpeedUp(data)

if __name__ == "__main__":
    main()
